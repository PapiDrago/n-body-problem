#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOLERANCE 1e-10

int main() {
    FILE *f1 = fopen("trajectory_serial.csv", "r");
    FILE *f2 = fopen("trajectory_parallel.csv", "r");

    if (!f1 || !f2) {
        printf("Error opening files.\n");
        return 1;
    }

    double val1, val2;
    int idx = 0, line = 1;
    int mismatch = 0;

    while (1) {
        int ret1 = fscanf(f1, "%lf", &val1);
        int ret2 = fscanf(f2, "%lf", &val2);

        if (ret1 != ret2) {
            printf("File length mismatch.\n");
            mismatch = 1;
            break;
        }

        if (ret1 == EOF) {
            break; // End of both files
        }

        double diff = fabs(val1 - val2);
        if (diff > TOLERANCE) {
            printf("Mismatch at line %d, element %d: %lf vs %lf (diff = %e)\n", line, idx, val1, val2, diff);
            mismatch = 1;
        }

        idx++;

        // Check for end of line in CSV (comma or newline)
        int c = fgetc(f1);
        if (c == '\n') {
            line++;
            idx = 0;
        }

        c = fgetc(f2);
        if (c == '\n') {
            // Already handled line++
        }
    }

    if (!mismatch) {
        printf("✅ Test passed: Files match within tolerance.\n");
    } else {
        printf("❌ Test failed: Differences found.\n");
    }

    fclose(f1);
    fclose(f2);
    return 0;
}

