#!/bin/sh

# This test generate N random sequence pairs of length L and an average error E using generate_dataset tool,
#   and pass them to quicked_harness.
# This test does not check for correctness, only for crashes.

# --- Cleanup ---

tempdir=$(mktemp -d)
cleanup() {
  rm -rf "$tempdir"
}

trap cleanup EXIT INT

# --- Input ---

L=$1
N=$2
E=$3

if [ -z "$N" ]; then
    N=10
fi

if [ -z "$L" ]; then
    L=10000
fi

if [ -z "$E" ]; then
    E=0.1
fi

# --- Test ---

echo "Generating $N random sequence pairs of length $L"
echo "$BIN_DIR"
"$BIN_DIR"/generate_dataset -l "$L" -n "$N" -e "$E" > "$tempdir/random_dataset.seq"

echo "Trimming sequences"
# Remove the first character of each line in random_dataset.seq
sed -i 's/^.//' "$tempdir/random_dataset.seq"

# For each sequence pair in random_dataset.seq, run ./quicked_harness seq1 seq2
# If quicked_harness crashes, abort the test with a non-zero exit code
echo "Running quicked_harness"
while read -r seq1 && read -r seq2; do
    # echo "Running quicked_harness with sequences:"
    # echo "> $seq1"
    # echo "< $seq2"
    # printf "\n"

    if ! "$BIN_DIR"/quicked_harness "$seq1" "$seq2"; then
        echo "quicked_harness crashed"
        exit 1
    fi
done < "$tempdir/random_dataset.seq"