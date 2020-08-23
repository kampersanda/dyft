#include <tinyformat/tinyformat.h>

static uint32_t HDT[256][256];

void compute_hdt(uint32_t bits) {
    const uint32_t N = 8 / bits;
    const uint32_t U = 1 << (N * bits);
    const uint32_t X = (1 << bits) - 1;

    for (uint32_t i = 0; i < 256; i++) {
        for (uint32_t j = 0; j < 256; j++) {
            HDT[i][j] = 0;
        }
    }

    for (uint32_t i = 0; i < U; i++) {
        for (uint32_t j = 0; j < U; j++) {
            uint32_t d = 0;
            for (uint32_t k = 0; k < N; k++) {
                const uint32_t a = (i >> (k * bits)) & X;
                const uint32_t b = (j >> (k * bits)) & X;
                d += (a != b) ? 1 : 0;
            }
            HDT[i][j] = d;
        }
    }
}

static uint8_t LUT[256][256];
static uint32_t DB[8][9];

void compute_lut(uint32_t bits) {
    const uint32_t N = 8 / bits;
    const uint32_t U = 1 << (N * bits);
    const uint32_t X = (1 << bits) - 1;

    for (uint32_t k = 0; k < 9; k++) {
        DB[bits - 1][k] = 0;
    }
    for (uint32_t j = 0; j < U; j++) {
        uint32_t d = 0;
        for (uint32_t k = 0; k < N; k++) {
            const uint32_t b = (j >> (k * bits)) & X;
            d += (0 != b) ? 1 : 0;
        }
        DB[bits - 1][d] += 1;
    }

    for (uint32_t i = 0; i < U; i++) {
        std::pair<uint32_t, uint32_t> row[256];  // c, d
        for (uint32_t j = 0; j < U; j++) {
            uint32_t d = 0;
            for (uint32_t k = 0; k < N; k++) {
                const uint32_t a = (i >> (k * bits)) & X;
                const uint32_t b = (j >> (k * bits)) & X;
                d += (a != b) ? 1 : 0;
            }
            row[j] = std::pair<uint32_t, uint32_t>{j, d};
        }
        std::sort(row, row + U, [](auto x, auto y) {
            if (x.second != y.second) {
                return x.second < y.second;
            }
            return x.first < y.first;
        });
        for (uint32_t j = 0; j < U; j++) {
            LUT[i][j] = row[j].first;
        }
        for (uint32_t j = U; j < 256; j++) {
            LUT[i][j] = 0;
        }
    }
}

int main(int argc, char** argv) {
#ifndef NDEBUG
    tfm::warnfln("The code is running in debug mode.");
#endif

    tfm::printfln("static constexpr uint8_t HD_TABLE[8][256][256] = {");
    for (uint32_t b = 1; b <= 8; b++) {
        compute_hdt(b);

        tfm::printfln("\t{ // b=%d", b);
        for (uint32_t i = 0; i < 256; i++) {
            tfm::printf("\t\t{");
            for (uint32_t j = 0; j < 256; j++) {
                tfm::printf("%d, ", HDT[i][j]);
            }
            tfm::printfln("},");
        }
        tfm::printfln("\t},");
    }
    tfm::printfln("};");

    tfm::printfln("static constexpr uint8_t LU_TABLE[8][256][256] = {");
    for (uint32_t b = 1; b <= 8; b++) {
        compute_lut(b);

        tfm::printfln("\t{ // b=%d", b);
        for (uint32_t i = 0; i < 256; i++) {
            tfm::printf("\t\t{");
            for (uint32_t j = 0; j < 256; j++) {
                tfm::printf("%d, ", LUT[i][j]);
            }
            tfm::printfln("},");
        }
        tfm::printfln("\t},");
    }
    tfm::printfln("};");

    tfm::printfln("static constexpr int BP_TABLE[8][10] = {");
    for (uint32_t b = 1; b <= 8; b++) {
        int n = 0;
        tfm::printf("\t{");
        for (uint32_t k = 0; k < 9; k++) {
            tfm::printf("%d, ", n);
            n += DB[b - 1][k];
        }
        tfm::printfln("%d}, // b=%d", n, b);
    }
    tfm::printfln("};");

    return 0;
}