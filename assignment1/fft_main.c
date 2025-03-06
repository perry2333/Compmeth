#include <stdio.h>
#include <stdlib.h>
#include "fft.c"  // Include the header file, not fft.c

int main() {
    int N_values[] = {64, 256, 1024, 4096};  // Corrected array declaration
    int testing_values[] = {0, 1, 2};
    int i, j;

    for (i = 0; i < 4; i++) {  // Fixed loop syntax
        for (j = 0; j < 3; j++) {
            char N_str[10], test_str[10];
            sprintf(N_str, "%d", N_values[i]);   // Convert integer to string
            sprintf(test_str, "%d", testing_values[j]);

            char *args[] = {"fft", N_str, test_str, NULL};  // Command-line arguments
            FFT_main(3, args);  // Call the fft function as if it were main()
        }
    }

    return 0;
}
