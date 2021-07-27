window.binty = (function() {
  const transform = {
    toCartesian([r, p, t]) {
      return [
          r*Math.sin(p)*Math.sin(t),
          r*Math.cos(p),
          r*Math.sin(p)*Math.cos(t),
      ];
    },
    toNoisedSpherical(coords, [nr, np, nt], radiusScaling = false) { // Note: not injecting many params in this file yet
      const [[radiusMin, radiusMax], [pMin, pMax], [tMin, tMax]] = getDimensionRanges(coords);
      const rNorm = radiusMax-radiusMin;

      // NOTE: scaling noise along the radius can benefit very low bitrates in partitioned mode,
      // but it isn't intuitive or attractive for geometric visual effects, so it's optional;
      // it should probably be invoked by default in partitioned compression and left out of the others
      return ([r, p, t], i, a) => {
          /*let rp =  p + (Math.random() * np * 2 - np) * r/rNorm;
          if (rp < 0) rp += Math.PI;
          if (rp > Math.PI) rp -= Math.PI;
          let rt =  t + (Math.random() * nt * 2 - nt) * r/rNorm;
          if (rt < -Math.PI) rt += 2*Math.PI;
          if (rt > Math.PI) rt -= 2*Math.PI;*/
          p = Math.max(Math.min(p, Math.PI - np*r/rNorm - Math.random()/20), 0 + np*r/rNorm + Math.random()/20);
          t = Math.max(Math.min(t, Math.PI - nt*r/rNorm - Math.random()/20), -Math.PI + nt*r/rNorm + Math.random()/20);

          return [
              r + (Math.random() * nr * 2 - nr) * (radiusScaling ? r/rNorm : 1),
              p + (Math.random() * np * 2 - np) * r/rNorm,
              t + (Math.random() * nt * 2 - nt) * r/rNorm,
          ];
      }
    },
    toNoisedCartesian(noise) {
      return Array.isArray(noise) ?
          (n, i, a) => n + Math.random()*noise[i%3]*2 - noise[i%3] :
          (n, i, a) => n + Math.random()*noise*2 - noise;
    },
    fromCustomTriplet([a, b, c]) {
      const [oa,ob,oc] = [a, b, c]
       .map((n, i) => [n,i])
       .sort((a, b) => a[0] - b[0])
       .map(n => n[1]);
      return n => [n[oa], n[ob], n[oc]];
    }
  };


  // Reads the final version of the file
  function fullRead(nums) { // TESTING ONLY! smartTruncate should be in HEADER!!!!!
    const TRIPLET_ORDER = [2,0,1]; // using constant naming convention; TODO: put this in header.

    // 1. decode 3bit and 5bit headers (4 bytes), of lengths 5 and 3 respectively. Remember to add 1 for them
    console.log('Starting Decompression.');
    const startTime = Date.now();
    const [mb1, mb2, mb3, colorBits, divCount] = fastMultiDecode(nums.slice(0,2),3, null, 5).map(n => n + 1);
    const multiBits = [mb1, mb2, mb3];
    const partitionLengths  = fastMultiDecode(nums.slice(2,4),5, null, 3).map(n => n + 1);

    // 2. decode fixed-length 32bits headers (6*4 noise + 6*4 for partitionLengths, plus 2*4 for color ranges = 56 bytes)
    const [
            nr, np, nt, nx, ny, nz,
            posLength, colorLength, // note: colorLength is de-alpha'd RGB length
            truncPosLength, truncColorLength,
            encodedPosLength, encodedColorLength,
            colorMin, colorMax,
    ] = unpackFloats(nums.slice(4, 4 + 56));
    const multiNoise = [nr, np, nt];
    const postNoise = [nx, ny, nz];
    const colorRange = [colorMin, colorMax];
    // 3. Use partitionLengths to figure out length of rest of the header
    const mangleSourcesLength = partitionLengths.reduce((acc, n) => acc + n, 0);
    const mangleSources = unpackFloats(nums.slice(4 + 56, 4 + 56 + 2*4*mangleSourcesLength));
    // 4. Reconstitute mangleMaps from last part of header
    const mangleMaps = buildMangleMapsFromFloats(mangleSources, partitionLengths, multiBits);
    // 5. Use divCount and multiBits to determine length of the positions array, run partitionDecode on it
    const encodedPositions = nums.slice(4 + 56 + 2*4*mangleSourcesLength, 4 + 56 + 2*4*mangleSourcesLength + encodedPosLength);
    const truncPositions = partitionDecode(encodedPositions, multiBits, partitionLengths[0], mangleMaps, truncPosLength);

    // 6. The rest of the array are the colors; run fastMultiDecode on it, using the color ranges decoded in #2
    const encodedColors = nums.slice(4 + 56 + 2*4*mangleSourcesLength + encodedPosLength);
    const colorDecoder = getMangleMaps(colorRange, colorBits, false, true)[1]; // isAlreadyRange true
    const truncColors = fastMultiDecode(encodedColors, colorBits, colorDecoder, truncColorLength);

    // 7. Detruncate and realpha
    //function smartDetruncate(nums, divCount, oLength, secondaries, secondaryUnitCount = 4) {

    const resultColors = realpha(detruncate(truncColors, colorLength)); // colorLength is rgb length; detruncate before realpha
    const fullPositions = detruncate(truncPositions, posLength);

    // 8. Use the noise data to add spherical noise, convert to cartesian, and add cartesian noise.
    const resultPositions = mapToCartesian(
        mapToNoisedSpherical(
            fastSplit(fullPositions, 3).map(transform.fromCustomTriplet(TRIPLET_ORDER)).flat(), // TODO: put TRIPLET_ORDER in header!
            multiNoise,
            true, // enable radius-scaling for radial noise; preferred for partition noise
        ),
    ).map(transform.toNoisedCartesian(postNoise));

    // 9. Convert to Float32Array and return.
    const result = [
        Float32Array.from(resultPositions),
        Float32Array.from(resultColors),
    ];
    console.log(`Decompression elapsed: ${Date.now() - startTime}ms`);
    return result;
  }

  function buildMangleMapsFromFloats(floats, partitionLengths, multiBits) {
    const output = [[],[],[]];
    let cursor = 0;
    for (let i = 0; i < partitionLengths.length; i++) {
        for (let j = 0; j < partitionLengths[i]; j++) {
            output[i][j] = getMangleMaps([floats[cursor++], floats[cursor++]], multiBits[i], false, true);
        }
    }
    return output;
  }

  // Combined from fastDecode and fastEulerDecode; seems just as fast for either
  function fastMultiDecode(nums, mb, ds, len, basis = 8) {
    let simple = typeof mb === 'number';
    if (!ds) ds = simple ? n => n : Array(3).fill(n => n);
    const length = nums.length;
    const bitArrayLength = length*basis;
    const bitArray = Array(bitArrayLength);
    for (let i = 0; i < length; i++) {
        const byte = (nums[i] >>> 0).toString(2).padStart(basis, '0');
        for (let j = 0; j < basis; j++) {
            bitArray[i*basis+j] = byte[j];
        }
      }
    output = [];
    const chunk = [];
    if (simple) {
      const sigArrayLength = mb*len;
      for (let i = 0; i <= sigArrayLength; i++) {
          if (i && i%mb === 0) output[i/mb-1] = ds(parseInt(chunk.join(''), 2));
          chunk[i%mb] = bitArray[i];
      }
      return output;
    }
    const bitSum = mb[0]+mb[1]+mb[2];
    // bitArrayLength minus remainder; assumes all triplets present; take ceiling in case missing
    const sigArrayLength = Math.ceil(bitSum*len/3);
    let outputLength = 0;
    let chunkCursor = 0;
    for (let i = 0; i <= sigArrayLength; i++) {
        const subIndex = i%bitSum;
        const a = (!subIndex || subIndex > mb[0] + mb[1]) ? 2 : (subIndex <= mb[0]) ? 0 : 1;
        const bits = mb[a];
        if (i && [mb[0], mb[0] + mb[1], 0].includes(subIndex)) {
            chunk.length = bits; // truncate chunk
            output[outputLength++] = ds[a](parseInt(chunk.join(''), 2));
            chunkCursor = 0;
        }
        chunk[chunkCursor++] = bitArray[i];
    }
    return output;
  }

  function mapPartitionSorted(nums, partitionCount = 2, cb = n => n) {
    const { ceil, floor, min } = Math;
    const p = partitionCount;
    let length = floor(nums.length/3);
    const orl = ceil(length/p);
    const opl = ceil(orl/p);
    const otl = ceil(opl/p);
    let outerIndex = 0;
    const output = Array(nums.length);
    for (let rp = 0; rp < p; rp++) {
        const rl = min(orl, length-rp*orl);
        for (let pp = 0; pp < p; pp++) {
            const pl = min(ceil(rl/p), rl-ceil(rl/p)*pp);
            for (tp = 0; tp < p; tp++) {
                const tl = Math.max(min(ceil(pl/p), pl-ceil(pl/p)*tp), 0);
                for (let i = 0; i < tl; i++) {
                    for (let dimensionIndex = 0; dimensionIndex < 3; dimensionIndex++) {
                        output[outerIndex] = cb([rp, pp, tp][dimensionIndex], nums[outerIndex], dimensionIndex, outerIndex);
                        outerIndex++;
                    }
                }

            }
        }
    }
    return output;
  }

  function realpha(rgbs) {
    const length = rgbs.length;
    const rgbas = [];
    let cursor = 0;
    for (let i = 0; i < length; i++) {
        rgbas[cursor++] = rgbs[i];
        if ((i+1)%3 === 0) rgbas[cursor++] = 1;
    }
    return rgbas;
  }

  function detruncate(divNums, oLength) {
    const output = Array(oLength);
    const divLength = divNums.length;
    let outputCursor = 0;
    for (let i = 0; i < oLength; i++) {
        output[outputCursor++] = divNums[i%divLength];
    }
    return output;
  }

  function mapToCartesian(coords) {
    return fastMapTriplets(coords, transform.toCartesian);
  }

  // Could maybe generalize to n-length chunks
  function fastMapTriplets(coords, mapper) {
    const remainder = coords.length%3;
    const length = coords.length;
    const triplet = [];
    const mapped = Array(coords.length - remainder); // truncate incomplete triplets
    let mappedLength = 0;
    for (let i = 0; i <= length; i++) {
        if (i && i%3 === 0) {
            const mappedTriplet = mapper(triplet);
            for (let j = 0; j < 3; j++) {
                mapped[mappedLength++] = mappedTriplet[j];
            }
        }
        triplet[i%3] = coords[i];
    }
    return mapped;
  }


  function fastSplit(nums, length) {
    const numsLength = nums.length;
    const outputLength = Math.ceil(numsLength/length);
    const output = Array(outputLength);
    const remainder = numsLength%length;
    let outputCursor = 0;
    let chunk, chunkCursor;
    for (let i = 0; i < numsLength; i++) {
        if (i%length === 0) {
            chunk = Array(length);
            chunkCursor = 0;
            output[outputCursor++] = chunk;
        }
        chunk[chunkCursor++] = nums[i];
    }
    if (remainder) {
     chunk.length = chunk.length - length + remainder; // include remainder, but remove overflow
    }
    return output;
  }

  function mapToNoisedSpherical(coords, multiNoise, radiusScaling = false) {
    return fastMapTriplets(coords, transform.toNoisedSpherical(coords, multiNoise, radiusScaling));
  }

  function getDimensionRanges(nums) {
    let l = nums.length;
    let max = [-Infinity, -Infinity, -Infinity];
    let min = [Infinity, Infinity, Infinity];
    while (l--) {
        let a = l%3;
        max[a] = nums[l] > max[a] ? nums[l] : max[a];
        min[a] = nums[l] < min[a] ? nums[l] : min[a];
    }
    return [[min[0], max[0]], [min[1], max[1]], [min[2], max[2]]];
  }


  function unpackFloats(ints) {
    return new Float32Array(Uint8Array.from(ints).buffer);
  }

  function getMangleMaps(nums, bits, includeSourceInfo = false, isAlreadyRange = false) {
    if(!nums.length) {
        console.warn('getMangleMaps called with empty array. Returning noop funcs.');
        return [() => {}, () => {}];
    }
    if (bits == undefined) {
        throw new Error('getMangleMaps requires a bit parameter.')
    }
    const [min, max] = isAlreadyRange ? nums : getRange(nums); // note: ranges are sometimes being sent in, making this redundant
    //console.log('num and range:', nums, min, max)
    if (max === min)  return maybeAddSourceInfo([n => max, n => max]);
    // make sure sub-unit ranges are visible after rounding:
    const vizFactor = _calcViz(max-min);
    //console.log(max-min, vizFactor)
    const bitMax = 2**bits-1;
    const compressMap = n => Math.round((n-min)/(max-min) * bitMax);
    const decompressMap = n => Math.round((n/bitMax * (max-min) + min)*vizFactor)/vizFactor;
    return maybeAddSourceInfo([compressMap, decompressMap]);

    function _calcViz(delta) {
         if (delta === -Infinity) {
             throw new Error('_calcViz called with Infinity delta; failing.');
         }
         let vizFactor = 1;
         while (delta * vizFactor < 100) {
             vizFactor*=10;
         }
         return vizFactor;
    }
    function maybeAddSourceInfo(result) {
        if (!includeSourceInfo) return result;
        return [...result, [min, max, bits]];
    }
  }

  function partitionDecode(nums, multiBits, partitionCount, partitionMangleMaps, length) {
    const firstPass = fastMultiDecode(
        nums,
        multiBits,
        [n => n, n => n, n => n], // passthrough decoder
        length,
    );
    const secondPass = mapPartitionSorted(
        firstPass,
        partitionCount,
        (encoderGroup, num, di) => partitionMangleMaps[di][encoderGroup][1](num),
    );
    return secondPass;
  }

  function loadBinary(filePath) {
    return fetch(filePath)
        .then(resp => resp.arrayBuffer())
        .then(buffer => new Uint8Array(buffer))
  }

  function applyData(positions, colors, posParent = _stars, colorParent = _starColors) {
    if (positions) {
        posParent.array = positions;
        posParent.needsUpdate = true;
    }
    if (colors) {
        colorParent.needsUpdate = true;
        colorParent.array = colors;
    }
  }
  return {
    load(filename) {
      return loadBinary(filename).then(fullRead);
    }
  }
})();