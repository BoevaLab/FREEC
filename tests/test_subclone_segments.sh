#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
build=$(mktemp -d)
trap 'rm -rf "$build"' EXIT
cxx=${CXX:-g++}
# Override CXXFLAGS to run without sanitizers on unsupported compilers.
read -r -a flags <<< "${CXXFLAGS:--O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer}"
sources=()
for source in "$repo"/src/*.cpp; do
    [[ $source == "$repo/src/main.cpp" ]] || sources+=("$source")
done
"$cxx" "${flags[@]}" -I"$repo/src" "$repo/tests/subclone_segments.cpp" "${sources[@]}" -lpthread -o "$build/subclone_segments"
"$build/subclone_segments" "$build/constant" constant
"$build/subclone_segments" "$build/segmented" segmented
