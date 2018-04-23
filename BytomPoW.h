#ifndef BYTOMPOW_H
#define BYTOMPOW_H

extern"C" {
    #include<cblas.h>  
}
#include "scrypt.h"
#include "sha3.h"
#include <iostream>
#include <assert.h>
#include <vector>
#include <stdint.h>

#define FNV(v1,v2) int32_t( ((v1)*FNV_PRIME) ^ (v2) )
const int FNV_PRIME = 0x01000193;

struct Mat256x256i8 {
    int8_t d[256][256];

    void toIdentityMatrix() {
        for(int i=0; i<256; i++) {
            for(int j=0; j<256; j++) {
                d[j][i]=0;
            }
        }
        for(int i=0; i<256; i++) {
            d[i][i]=1;
        }
    }

    void copyFrom(const Mat256x256i8& other) {
        for(int i=0; i<256; i++) {
            for(int j=0; j<256; j++) {
                this->d[j][i]=other.d[j][i];
            }
        }
    }

    Mat256x256i8() {
        this->toIdentityMatrix();
    }

    Mat256x256i8(const Mat256x256i8& other) {
        this->copyFrom(other);
    }

    void copyFrom_helper(LTCMemory& ltcMem, int offset) {
        for(int i=0; i<256; i++) {
            const Words32& lo=ltcMem.get(i*4+offset);
            const Words32& hi=ltcMem.get(i*4+2+offset);
            for(int j=0; j<64; j++) {
                uint32_t i32=j>=32?hi.get(j-32):lo.get(j);
                d[j*4+0][i]=(i32>>0)&0xFF;
                d[j*4+1][i]=(i32>>8)&0xFF;
                d[j*4+2][i]=(i32>>16)&0xFF;
                d[j*4+3][i]=(i32>>24)&0xFF;
            }
        }
    }

    void copyFromEven(LTCMemory& ltcMem) {
        copyFrom_helper(ltcMem, 0);
    }

    void copyFromOdd(LTCMemory& ltcMem) {
        copyFrom_helper(ltcMem, 1);
    }

    void mul(const Mat256x256i8& a, const Mat256x256i8& b) {
        double *ma = (double *)calloc(256*256, sizeof(double));
        double *mb = (double *)calloc(256*256, sizeof(double));
        double *mc = (double *)calloc(256*256, sizeof(double));

        for (int i = 0; i < 256; i++) {
            for (int j = 0; j < 256; j++) {
                ma[i*256+j] = (double)(a.d[i][j]);
                mb[i*256+j] = (double)(b.d[i][j]);
            }
        }

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 256, 256, 256, 1, ma, 256, mb, 256, 0, mc, 256);

        for (int i = 0; i < 256; i++) {
            for (int j = 0; j < 256; j++) {
                int tmp = (int64_t)(mc[i*256+j]);
                this->d[i][j]=((tmp&0xFF)+ ((tmp>>8)&0xFF))&0xFF;
            }
        }

        free(ma);
        free(mb);
        free(mc);
    }

    void add(Mat256x256i8& a, Mat256x256i8& b) {
        for(int i=0; i<256; i++) {
            for(int j=0; j<256; j++) {
                int tmp=int(a.d[i][j])+int(b.d[i][j]);
                this->d[i][j]=(tmp&0xFF);
            }
        }
    }
};

struct Arr256x64i32 {
    uint32_t d[256][64];

    uint8_t* d0RawPtr() {
        return (uint8_t*)(d[0]);
    }

    Arr256x64i32(const Mat256x256i8& mat) {
        for(int j=0; j<256; j++) {
            for(int i=0; i<64; i++) {
                d[j][i] = ((uint32_t(uint8_t(mat.d[j][i  + 192]))) << 24) |
                          ((uint32_t(uint8_t(mat.d[j][i + 128]))) << 16) |
                          ((uint32_t(uint8_t(mat.d[j][i  + 64]))) << 8) |
                          ((uint32_t(uint8_t(mat.d[j][i ]))) << 0);
            }
        }
    }

    void reduceFNV() {
        for(int k=256; k>1; k=k/2) {
            for(int j=0; j<k/2; j++) {
                for(int i=0; i<64; i++) {
                    d[j][i] = FNV(d[j][i], d[j + k / 2][i]);
                }
            }
        }
    }
};

struct BytomMatList {
    std::vector<Mat256x256i8*> matVec;

    Mat256x256i8 at(int i) {
        return *(matVec[i]);
    }

    BytomMatList() {
        for(int i=0; i<256; i++) {
            Mat256x256i8* ptr = new Mat256x256i8;
            assert(ptr!=NULL);
            matVec.push_back(ptr);
        }
    }

    ~BytomMatList() {
        for(int i=0; i<256; i++) {
            delete matVec[i];
        }
    }

    void init(const Words32& X_in) {
        Words32 X = X_in;
        LTCMemory ltcMem;
        for(int i=0; i<128; i++) {
            ltcMem.scrypt(X);
            matVec[2*i]->copyFromEven(ltcMem);
            matVec[2*i+1]->copyFromOdd(ltcMem);
        }
    }
};

extern BytomMatList* matList_int8;

static inline void iter_mineBytom(
                        const uint8_t *fixedMessage,
                        uint32_t len,
                        uint8_t result[32]) {
    Mat256x256i8 *res = new Mat256x256i8[4];
    Mat256x256i8 *mat=new Mat256x256i8;
    sha3_ctx *ctx = (sha3_ctx*)calloc(1, sizeof(*ctx));

    clock_t start, end;
    start = clock();
    for(int k=0; k<4; k++) { // The k-loop
        uint8_t sequence[128];
        rhash_sha3_256_init(ctx);
        rhash_sha3_update(ctx, fixedMessage+(len*k/4), len/4);//分四轮消耗掉fixedMessage
        rhash_sha3_final(ctx, sequence);
        Mat256x256i8 *tmp=new Mat256x256i8;
        tmp->toIdentityMatrix();
        for(int j=0; j<2; j++) {
            for(int i=0; i<32; i++) {
                mat->mul(*tmp, matList_int8->at(sequence[i]));
                tmp->copyFrom(*mat);
            }
        }
        res[k].copyFrom(*mat);
        delete tmp;
    }

    mat->add(res[0], res[1]);
    mat->add(*mat, res[2]);
    mat->add(*mat, res[3]);

    end = clock();
    std::cout << "\tTime for getting MulMatix: "
              <<(double)(end - start) / CLOCKS_PER_SEC << "s"
              << std::endl;

    Arr256x64i32 arr(*mat);
    arr.reduceFNV();
    rhash_sha3_256_init(ctx);
    rhash_sha3_update(ctx, arr.d0RawPtr(), 256);
    rhash_sha3_final(ctx, result);
    delete mat;
    delete[] res;
    free(ctx);
}

static inline void incrNonce(uint8_t nonce[8]) {
    for(int i=0; i<8; i++) {
        if(nonce[i]!=255) {
            nonce[i]++;
            break;
        } else {
            nonce[i]=0;
        }
    }
}

static inline int countLeadingZero(uint8_t result[32]) {
    int count=0;
    for (int i=31; i>=0; i--) { // NOTE: reverse
        if (result[i] < 1) {
            count+=8;
        } else if (result[i]<2)  {
            count+=7;
            break;
        } else if (result[i]<4)  {
            count+=6;
            break;
        } else if (result[i]<8)  {
            count+=5;
            break;
        } else if (result[i]<16) {
            count+=4;
            break;
        } else if (result[i]<32) {
            count+=3;
            break;
        } else if (result[i]<64) {
            count+=2;
            break;
        } else if (result[i]<128) {
            count+=1;
            break;
        }
    }
    return count;
}

#endif

