#!/usr/bin/env sh
set -vex

db_url="$1"
output="$2"
md5_hash="$3"

wget "$db_url" -O "$output"

md5=$(md5sum "$output" | awk '{ print $1 }')


if [ "$md5" = "$md5_hash" ]; then
    # The MD5 sum matched
    tar xzf "$output" -C $(dirname $output)
    rm "$output"
else
    # The MD5 sum didn't match
    echo "Centrifuge database download's md5sum did not match expected."
    echo "Got $md5"
    echo "Expected $md5_hash"
    exit 1
fi

exit 0
