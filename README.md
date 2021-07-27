### Binty Rust

#### Summary
This is a toy project for a naive integer compressor I'm using [here](https://shreet.surge.sh) to get a galaxy model down from a binary GLTF of 1.5mb to 30k and 170k (looks fine). In the linked example, it'll big-bang the 30k galaxy until the 170k one loads.

Unfortunately, the decompression takes about 150ms, and I wanted to see how fast it is in Rust/WASM.

#### Description

The gist of how the compressor works:

1. It assumes that we're just dealing with a big Float32 array of cartesian coordinates (as well as an array of RGBA values corresponding to those coordinates, but it's dealt with separately).
2. The compressor includes a lot of parameters for truncating and re-replicating parts of the array, whether to compress in cartesian or spherical coordinates, how many bits to use per dimension, and how many partitions to divide the coordinate space into.
3. Each partition gets 6 Float32s to represent its bounds.
4. The other coordinates in a partition are integer compressed and bitpacked according to the partition bounds (ie. 3bits would give a point 8 possible positions per dimension, oriented with respect to the partition bounds upon decompression).
5. The array is sorted in the following fashion, with the loader simply mapping partition bounds onto each coordinate using modulo (just the top half):

<p align="center">
    <img src="https://user-images.githubusercontent.com/2810775/127161534-c47da5e8-82c3-4de5-a351-15beb9570c5c.png"/>
</p>

6. The bottom half of the above shows interlacing of partitions for progressive loading. It's not currently included in the format, as I haven't quite finished it.
