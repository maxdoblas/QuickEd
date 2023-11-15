#include <quicked.h>
#include <stdio.h>

int main() {
    quicked_aligner_t aligner;
    quicked_params_t params = quicked_default_params();

    params.algo = BANDED;
    params.bandwidth = 50;
    params.only_score = true;

    quicked_new(&aligner, params);

    quicked_align(&aligner, "ACGT", 4, "ACTT", 4);

    printf("Score: %d\n", aligner.score);
    printf("CIGAR <Expecting NULL>: %s\n", aligner.cigar);

    quicked_free(&aligner);

    return 0;
}