

/*================================================= CODE STOLEN FROM MY OTHER PROJECT ======================*/

/* Deterministic 2d random */
function det_random_2d(seed, layers, layer, x, y) {

    // merge coords with layers
    if (x % 2 == y % 2) {
        x = x * layers;
        y = y * layers;
        x = x + layer;
    }
    else {
        x = x * layers;
        y = y * layers;
        y = y + layer;
    }
    
    //-------------------------------------------------------------------------
    
    // compute
    let val = 0;

    // temporary vars for linear transformations
    let nx, ny;

    x += seed;
    
    // xor shuffle
    y ^= x;
    nx = x - y;
    ny = x + y - seed;
    x = nx;
    y = ny;
    x ^= y;

    // radius shuffle
    y ^= ((x * x) & -Math.abs(x)) * Math.sign(x);
    x ^= ((y * y) & -Math.abs(y)) * Math.sign(y);
    y ^= ((x * x) & -Math.abs(x)) * Math.sign(x);
    x ^= ((y * y) & -Math.abs(y)) * Math.sign(y);
    
    val = x ^ y ^ seed;
   
    //-----------------------------------------------------------
    
    // and return
    if (Math.floor(val / 997) % 2 == 0)
        val = ((997 + (val % 997)) % 997) / 997;
    else
        val = ((997 + (-val % 997)) % 997) / 997;
       
    return val;
}

//==========================================================================================

const HALF_PI = Math.PI / 2;
const TWO_PI = Math.PI * 2;

function smoothf_norm(val) {
    return Math.sin(val * HALF_PI);
}

// smooth
function smoothf(val) {

    let fixval = val * 2 - 1;

    let smoothval = fixval;

    smoothval = smoothf_norm(smoothval);

    let backval = (smoothval + 1) / 2;

    return backval;
};

let pnoise_cache = [];
let pnoise_cache_entries = 2000;

for (let k = 0; k < pnoise_cache_entries; ++k) {
    pnoise_cache.push({ ci : 0, cj : 0, cmseed : 0, crsize : 0, value : 0});
}

function pnoise(i, j, mseed, rsize) {

    // check cache
    let hash = Math.abs(((rsize << 8) ^ (mseed << 4) ^ (i << 2) ^ (j) ^ (-1))) % pnoise_cache_entries;

    let def = pnoise_cache[hash];

    if (def.crsize === rsize && def.cmseed === mseed && def.ci === i && def.cj === j) {
        return def.value;
    }
    
    rsize = Math.floor(rsize);
    
    if (rsize === 0)
        return 0.5;
    
    i = Math.floor(i);
    j = Math.floor(j);

    let idrsize = Math.floor(i / rsize);
    let jdrsize = Math.floor(j / rsize);

    let tl_x, tl_y, bl_x, bl_y, tr_x, tr_y, br_x, br_y;

    let tl_r = det_random_2d(mseed, 1, 0, idrsize, jdrsize);
    let tr_r = det_random_2d(mseed, 1, 0, idrsize + 1, jdrsize);
    let bl_r = det_random_2d(mseed, 1, 0, idrsize, jdrsize + 1);
    let br_r = det_random_2d(mseed, 1, 0, idrsize + 1, jdrsize + 1);
            
    tl_x = Math.sin(TWO_PI * tl_r);
    tr_x = Math.sin(TWO_PI * tr_r);
    bl_x = Math.sin(TWO_PI * bl_r);
    br_x = Math.sin(TWO_PI * br_r);

    tl_y = Math.cos(TWO_PI * tl_r);
    tr_y = Math.cos(TWO_PI * tr_r);
    bl_y = Math.cos(TWO_PI * bl_r);
    br_y = Math.cos(TWO_PI * br_r);
    

    //---------------------------------------------------------------
    let ceili = Math.ceil(i / rsize);
    if (ceili == idrsize)
        ceili++;

    let ceilj = Math.ceil(j / rsize);
    if (ceilj == jdrsize)
        ceilj++;

    let idr_r_size = idrsize * rsize;
    let jdr_r_size = jdrsize * rsize;

    let dx = (i - idr_r_size) / (ceili * rsize - idr_r_size);
    let dy = (j - jdr_r_size) / (ceilj * rsize - jdr_r_size);

    let dxm = dx - 1;
    let dym = dy - 1;

    let tl = (dx * tl_x) + (dy * tl_y);
    let tr = (dxm * tr_x) + (dy * tr_y);
    let bl = (dx * bl_x) + (dym * bl_y);
    let br = (dxm * br_x) + (dym * br_y);

    dx = smoothf(dx);
    dy = smoothf(dy);

    tl = smoothf_norm(tl);
    tr = smoothf_norm(tr);
    bl = smoothf_norm(bl);
    br = smoothf_norm(br);

    dxm = 1 - dx;
    dym = 1 - dy;

    let power = (tl * dxm + tr * dx) * dym + (bl * dxm + br * dx) * dy;

    let val = (power + 1) / 2;

    def.ci = i;
    def.cj = j;
    def.cmseed = mseed;
    def.crsize = rsize;
    def.value = val;
    
    return val;
}

/*================================================= UTIL CODE ==============================================*/

/* Appends source contents to the end of target contents (in place!) */
function flush_array(target, source) {
    for (let i = 0; i < source.length; ++i) {
        target.push(source[i]);
    }
    return target;
}

/* Replicates array a few times */
function replicate_array(source, repeat) {
    let target = [];
    while (repeat !== 0) {
        flush_array(target, source);
        repeat--;
    }
    return target;
}

/* Utils for 2D UV manipulation */
class Math2dForUV {};

Math2dForUV.TOP_LEFT = "tl";
Math2dForUV.TOP_RIGHT = "tr";
Math2dForUV.BOTTOM_LEFT = "bl";
Math2dForUV.BOTTOM_RIGHT = "br";

/* normalize_coords is: TOP_LEFT, TOP_RIGHT, BOTTOM_RIGHT, BOTTOM_LEFT */
Math2dForUV.selectUV = function(selectors, normalizedCoords) {

    let results = [];

    for (let sel of selectors) {
        if (sel === Math2dForUV.TOP_LEFT) {
            results.push(normalizedCoords[0]);
            results.push(normalizedCoords[1]);
        }
        if (sel === Math2dForUV.TOP_RIGHT) {
            results.push(normalizedCoords[2]);
            results.push(normalizedCoords[3]);
        }
        if (sel === Math2dForUV.BOTTOM_RIGHT) {
            results.push(normalizedCoords[4]);
            results.push(normalizedCoords[5]);
        }
        if (sel === Math2dForUV.BOTTOM_LEFT) {
            results.push(normalizedCoords[6]);
            results.push(normalizedCoords[7]);
        }
    }
                
    return results;
}

/* Adds added array to a source, iterating over added array as a loop (in place!) */
function cyclic_add_array(source, added) {
    for (let i = 0; i < source.length; ++i) {
        let ai = i % added.length;
        source[i] += added[ai];
    }
    return source;
}

/* More than just Math class */
class BeyondMath {};

/* Returns one with "one" probability, otherwise returns zero */
BeyondMath.indicatorRandom = function (one) {
    if (Math.random() < one) return 1;
    return 0;
}

/*================================================= ENGINE-INDEPENDENT CLASSES ==============================================*/


class BlockUvMap {
    
    /* l r u d f b */
    project(idBlock, side) {

        if (idBlock === 3) {
            
            if (side === 'u') {
                return [0.5, 0, 1, 0, 1, 0.5, 0.5, 0.5];    // grass
            }

            if (side === 'd') {
                return [0.5, 0.5, 1, 0.5, 1, 1, 0.5, 1];    // dirt
            }

            return [0, 0.5, 0.5, 0.5, 0.5, 1, 0, 1];    // dirt grass
        }

        if (idBlock === 2) {
            return [0.5, 0.5, 1, 0.5, 1, 1, 0.5, 1];  // dirt
        }

        if (idBlock === 1) {
            return [0, 0, 0.5, 0, 0.5, 0.5, 0, 0.5]; // stone
        }

        return [0, 0, 1, 0, 1, 1, 0, 1];    // default - all atlas
    }
}

/* Represents a class for making custom terrain from squares */
class SquareGridBuilder {
    
    constructor() {
        this.positions =  [];
        this.indices = [];
        this.uvs = [];
        this.normals = [];
        this.uvMap = new BlockUvMap();
        
        
        this.straightIndices = [0, 2, 1, 1, 2, 3];
        this.reverseIndices = [1, 2, 0, 3, 2, 1]; // just like straight, but the alternative side
    }

    insertPiece(new_positions, new_indices, new_uvs, new_normals) {
        flush_array(this.positions, new_positions);
        flush_array(this.indices, new_indices);
        flush_array(this.uvs, new_uvs);
        flush_array(this.normals, new_normals);
    }
    
    insertSquareGridBuilder(otherSgb) {
        
        let all_indices = this.positions.length / 3;
        
        flush_array(this.positions, otherSgb.positions);
        flush_array(this.indices, cyclic_add_array([...otherSgb.indices], [all_indices]));
        flush_array(this.uvs, otherSgb.uvs);
        flush_array(this.normals, otherSgb.normals);
    }

    insertSide(xyz, vertex_loop, index_loop, normal_sample, uvs) {
        
        let all_indices = this.positions.length / 3;
        
        this.insertPiece(
            cyclic_add_array(vertex_loop, xyz)
            , cyclic_add_array(index_loop, [all_indices])                    
            , uvs
            , replicate_array(normal_sample, 4)
        );
    }

    insertLeft(x, y, z, idBlock) {
        
        this.insertSide(
            [x, y, z]
            , [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1]
            , [...this.straightIndices]
            , [-1, 0, 0]
            , Math2dForUV.selectUV(["tr", "tl", "br", "bl"], this.uvMap.project(idBlock, "l"))
        );
    }

    insertRight(x, y, z, idBlock) {
        
        this.insertSide(
            [x + 1, y, z]
            , [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1]
            , [...this.reverseIndices]
            , [1, 0, 0]
            , Math2dForUV.selectUV(["tl", "tr", "bl", "br"], this.uvMap.project(idBlock, "r"))
        );
    }

    insertBack(x, y, z, idBlock) {
        
        this.insertSide(
            [x, y, z]
            , [0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0]
            , [...this.straightIndices]
            , [0, 0, -1]
            , Math2dForUV.selectUV(["tl", "bl", "tr", "br"], this.uvMap.project(idBlock, "b"))
        );
    }

    insertFront(x, y, z, idBlock) {
        
        this.insertSide(
            [x, y, z + 1]
            , [0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0]
            , [...this.reverseIndices]
            , [0, 0, 1]
            //, [1, 0, 1, 1, 0, 0, 0, 1]
            , Math2dForUV.selectUV(["tr", "br", "tl", "bl"], this.uvMap.project(idBlock, "f"))
        );
    }

    insertTop(x, y, z, idBlock) {
        
        this.insertSide(
            [x, y + 1, z]
            , [0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1]
            , [...this.straightIndices]
            , [0, 1, 0]
            , Math2dForUV.selectUV(["tl", "bl", "tr", "br"], this.uvMap.project(idBlock, "u"))
        );
    }

    insertBottom(x, y, z, idBlock) {
        
        this.insertSide(
            [x, y, z]
            , [0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1]
            , [...this.reverseIndices]
            , [0, -1, 0]
            , Math2dForUV.selectUV(["tr", "br", "tl", "bl"], this.uvMap.project(idBlock, "d"))
        );
    }
}

class Array3D {
    
    constructor(sx, sy, sz) {
        this.sx = sx;
        this.sy = sy;
        this.sz = sz;
        this.data = [];
    }

    hasCoords(x, y, z) {
        return (x >= 0 && x < this.sx) && (y >= 0 && y < this.sy) && (z >= 0 && z < this.sz);
    }

    toLinear(x, y, z) {
        return z + this.sz * (y + x * this.sy);
    }

    setValue(x, y, z, value) {
        this.data[this.toLinear(x, y, z)] = value;
    }

    getValue(x, y, z) {
        return this.data[this.toLinear(x, y, z)];
    }

    getValueEx(x, y, z, otherwise) {
        if (!this.hasCoords(x, y, z))
            return otherwise;
        return this.data[this.toLinear(x, y, z)];
    }

    iterate(iter) {
        for (let i = 0; i < this.sx; ++i) {
            for (let j = 0; j < this.sy; ++j) {
                for (let k = 0; k < this.sz; ++k) {
                    let new_val = iter(i, j, k, this.getValue(i, j, k));
                    if (new_val !== undefined)
                        this.setValue(i, j, k, new_val);
                }
            }
        }
    }

    iterateOuterSide(iter) {

        let lowEnd = 0;
        let highEnd = this.sz - 1;
        
        for (let i = 0; i < this.sx; ++i) {

            let iEnd = i === lowEnd || i === highEnd;

            for (let j = 0; j < this.sy; ++j) {     

                let jEnd = j === lowEnd || j === highEnd;
                
                if (iEnd || jEnd) {
                    for (let k = 0; k < this.sz; ++k) {
                        let new_val = iter(i, j, k, this.getValue(i, j, k));
                        if (new_val !== undefined)
                            this.setValue(i, j, k, new_val);
                    }
                }
                else {
                    let new_val = iter(i, j, lowEnd, this.getValue(i, j, lowEnd));
                    if (new_val !== undefined)
                        this.setValue(i, j, lowEnd, new_val);

                    new_val = iter(i, j, highEnd, this.getValue(i, j, highEnd));
                    if (new_val !== undefined)
                        this.setValue(i, j, highEnd, new_val);
                }
            }
        }
        
        //--------                     
    }
    
}

const CHUNK_SIZE = 16;

/* Represents a single terrain chunk 16 * 16 * 16 size */
class TerrainChunk {

    /* TODO: make sparse chunks - like - all air */
    
    /* init chunks own coodrs */
    constructor(cx, cy, cz) {
        this.cx = cx;
        this.cy = cy;
        this.cz = cz;
        this.chunkData = new Array3D(CHUNK_SIZE, CHUNK_SIZE, CHUNK_SIZE);   // TODO: make loadable
        this.visible = true;
        this.visibleDirty = false;  // if any attempt to update it while invicible encountered

        this.resetDeep();
    }

    resetDeep() {
        this.innerSgb = null;
        this.outerSgb = null;
    }
    
    /* TODO: make chunk of central and peripheral sides, to optimize recalculations */
    reset() {
        //this.localSgb = null;

        //this.innerSgb = null;
        this.outerSgb = null;
    }
    
    drawToSGB(sgb, leftNgh, rightNgh, topNgh, bottomNgh, frontNgh, backNgh) {
        
        if (this.outerSgb) {
            sgb.insertSquareGridBuilder(this.innerSgb);
            sgb.insertSquareGridBuilder(this.outerSgb);
            return;
        }

        let noInnerSide = !this.innerSgb;

        if (noInnerSide)
            this.innerSgb = new SquareGridBuilder();
        
        this.outerSgb = new SquareGridBuilder();

        let outerSgb = this.outerSgb;
        let innerSgb = this.innerSgb;

        let offsetX = CHUNK_SIZE * this.cx;
        let offsetY = CHUNK_SIZE * this.cy;
        let offsetZ = CHUNK_SIZE * this.cz;

        let chunkData = this.chunkData;

        let iteratingFunction = (noInnerSide ? ((iter) => chunkData.iterate(iter)) : ((iter) => chunkData.iterateOuterSide(iter)));
                 
        iteratingFunction(
            (i, j, k, idBlock) => {

                let isOuter = (i === 0 || j === 0 || k === 0 || i === CHUNK_SIZE - 1 || j === CHUNK_SIZE - 1 || k === CHUNK_SIZE - 1);

                let localSgb = isOuter ? outerSgb : innerSgb;

                let airBlockCode = 0;

                if (idBlock !== airBlockCode) {
                    
                    /* TODO: replace === airBlockCode with a newly created isTransparent() function */
                    if (chunkData.getValueEx(i - 1, j, k, TerrainChunk.testForValueAt(leftNgh, CHUNK_SIZE - 1, j, k, airBlockCode)) === airBlockCode)
                        localSgb.insertLeft(i + offsetX, j + offsetY, k + offsetZ, idBlock);
                    
                    if (chunkData.getValueEx(i + 1, j, k, TerrainChunk.testForValueAt(rightNgh, 0, j, k, airBlockCode)) === airBlockCode)
                        localSgb.insertRight(i + offsetX, j + offsetY, k + offsetZ, idBlock);
                    
                    if (chunkData.getValueEx(i, j + 1, k, TerrainChunk.testForValueAt(topNgh, i, 0, k, airBlockCode)) === airBlockCode)
                        localSgb.insertTop(i + offsetX, j + offsetY, k + offsetZ, idBlock);
                    
                    if (chunkData.getValueEx(i, j - 1, k, TerrainChunk.testForValueAt(bottomNgh, i, CHUNK_SIZE - 1, k, airBlockCode)) === airBlockCode)
                        localSgb.insertBottom(i + offsetX, j + offsetY, k + offsetZ, idBlock);
                    
                    if (chunkData.getValueEx(i, j, k + 1, TerrainChunk.testForValueAt(frontNgh, i, j, 0, airBlockCode)) === airBlockCode)
                        localSgb.insertFront(i + offsetX, j + offsetY, k + offsetZ, idBlock);
                    
                    if (chunkData.getValueEx(i, j, k - 1, TerrainChunk.testForValueAt(backNgh, i, j, CHUNK_SIZE - 1, airBlockCode)) === airBlockCode)
                        localSgb.insertBack(i + offsetX, j + offsetY, k + offsetZ, idBlock);
                    
                }
            }
        );
        
        sgb.insertSquareGridBuilder(this.innerSgb);
        sgb.insertSquareGridBuilder(this.outerSgb);
    }
}

TerrainChunk.testForValueAt = function (tc, x, y, z, otherwise) {
    if (tc === null || tc === undefined)
        return otherwise;            
    return tc.chunkData.getValueEx(x, y, z, otherwise);
}

/* Class responsible for loading chunks from particular coords */
class TerrainChunkProvider {
    
    constructor() {}

    obtainChunkFrom(cx, cy, cz) {

        let chunk = new TerrainChunk(cx, cy, cz);

        chunk.chunkData.iterate(
            (x, y, z, value) => {

                //f ((cx ^ cy ^ cz) % 8 !== 0) return 0;
                let gx = CHUNK_SIZE * cx + x;
                let gy = CHUNK_SIZE * cy + y;
                let gz = CHUNK_SIZE * cz + z;

                //return 1;
                let height = pnoise(gx, gz, 1337666, 20);

                if (height > 25) return 0;
                if (height < -25) return 1;

                let maxHeight = Math.floor(50 * height) - 25;

                if (gy > maxHeight) return 0;
                if (gy === maxHeight) return 3;
                if (gy > maxHeight - 3) return 2;
                
                return 1;
            }
        );
        
        return chunk;
    }
}

/* Class responsible for chuck loading and unloading, and their cumulative rendering with the neighbour consideration */
class TerrainGrid {
    
    constructor(renderer) {
        this.chunkMap = new Object();
        this.changed = false;
        this.renderer = renderer;
    }

    dropChunk(x, y, z) {

        let cm = this.chunkMap;
        let droppedChunk = cm[[x, y, z]];

        delete cm[[x, y, z]];

        if (droppedChunk !== null && droppedChunk !== undefined)
            this.updateChunkNeighboursAt(x, y, z);

        this.renderer.deleteChunk(x, y, z);
    }
    
    uploadChunk(chunk) {
        //console.log("TERRAIN: upload " + [chunk.cx, chunk.cy, chunk.cz]);
        
        let cm = this.chunkMap;
        chunk.reset();
        cm[[chunk.cx, chunk.cy, chunk.cz]] = chunk;

        this.updateChunkNeighboursAt(chunk.cx, chunk.cy, chunk.cz);
        this.updateChunkAt(chunk.cx, chunk.cy, chunk.cz);
    }

    updateVisibility(x, y, z, visible) {
        let ch = this.getChunk(x, y, z);
        
        /* invisible chunks might have become dirty for long no-updates, let's refresh them */
        if (ch) {

            if (!ch.visible && visible && ch.visibleDirty) {
                ch.visible = visible;
                this.updateChunkAt(x, y, z);
            }
            else {
                ch.visible = visible;
            }
            
            this.renderer.updateVisibility(x, y, z, visible);
        }
    }

    getChunk(x, y, z) {
        let cm = this.chunkMap;
        return cm[[x, y, z]];
    }

    hasChunk(x, y, z) {
        let gc = this.getChunk(x, y, z);
        return !(gc === undefined || gc === null); // TODO: check this - can I use just "!!gc" ?
    }

    updateChunkAt(x, y, z) {
        let cn = this.getChunk(x, y, z);

        if (cn !== null && cn !== undefined) {

            /* If this chunk is invicible, why should we ever update it's structure? It's computationally expensive */
            if (cn.visible === false) {
                cn.visibleDirty = true;
                return;
            }

            cn.reset();

            let cn_sgb = new SquareGridBuilder();
                
            cn.drawToSGB(
                cn_sgb
                , this.getChunk(x - 1, y, z)
                , this.getChunk(x + 1, y, z)
                , this.getChunk(x, y + 1, z)
                , this.getChunk(x, y - 1, z)
                , this.getChunk(x, y, z + 1)
                , this.getChunk(x, y, z - 1)                        
            );

            this.renderer.updateChunk(x, y, z, cn_sgb);
        }
    }

    updateChunkNeighboursAt(x, y, z) {

        this.updateChunkAt(x - 1, y, z);
        this.updateChunkAt(x + 1, y, z);
        this.updateChunkAt(x, y - 1, z);
        this.updateChunkAt(x, y + 1, z);
        this.updateChunkAt(x, y, z - 1);
        this.updateChunkAt(x, y, z + 1);

        this.changed = true;
    }

}

//const CHUNKS_VISIBLE_RADIUS_NEAR = 2;
const CHUNKS_VISIBLE_RADIUS = 8;    // marked for addition
const CHUNKS_VISIBLE_HEIGHT_RADIUS = 4; // marked for addition height only

const CHUNKS_VISIBLE_RADIUS_FAR = 10;   // marked for deletion
const CHUNKS_VISIBLE_RADIUS_PSEUDO_FAR = 8;    // not visible even if not deleted yet

const DEFAULT_BATCH_LEVEL = 0;

/* TODO: I need a class that would manually construct TerrainGrid of chunks and update new chunks in it as the new data approaches */
class TerrainGridManipulator {

    constructor(terrain, chunkProvider) {
        this.terrain = terrain;
        this.chunkProvider = chunkProvider;

        this.pendingAsyncTasks = new Object();
        
        this.asyncRunning = false;
        this.dominatingChunkActionPriority = 1;

        this.cx = this.cy = this.cz = undefined;
        this.smallCx = this.smallCy = this.smallCz = undefined;
    }
    
    generate() {
        let terrain = this.terrain;
        let chunkProvider = this.chunkProvider;
    }

    cancelTaskAt(x, y, z) {
        delete this.pendingAsyncTasks[[x, y, z]];
    }

    orderTaskAt(x, y, z, task) {
        this.pendingAsyncTasks[[x, y, z]] = [x, y, z, task];
        this.kickstartAsyncChain();
    }

    getTaskAt(x, y, z) {
        return this.pendingAsyncTasks[[x, y, z]]?.[3];
    }

    askAsyncUploadAt(x, y, z) {

        if (!this.terrain.hasChunk(x, y, z)) {
            this.orderTaskAt(x, y, z, 1);
        }
        else {
            if (this.getTaskAt(x, y, z) === 0) {
                this.cancelTaskAt(x, y, z);
            }
        }
    }

    askAsyncDeleteAt(x, y, z) {

        if (this.terrain.hasChunk(x, y, z)) {
            this.orderTaskAt(x, y, z, 0);
        }
        else {
            if (this.getTaskAt(x, y, z) === 1) {
                this.cancelTaskAt(x, y, z);
            }
        }
    }

    refreshTaskPriority() {

        this.addTasks = [];
        this.removeTasks = [];

        for (const task of Object.entries(this.pendingAsyncTasks)) {

            let taskType = task[1][3];

            if (taskType === 0)
                this.removeTasks.push(task);

            if (taskType === 1)
                this.addTasks.push(task);
        }

        let that = this;

        this.addTasks = this.addTasks.sort(function (left, right) {
            let leftDef = left[1];
            let rightDef = right[1];

            let leftDistance = 2 * Math.abs(leftDef[0] - that.cx) + /* Math.abs(leftDef[1] - that.cy) + */ 2 * Math.abs(leftDef[2] - that.cz);
            let rightDistance = 2 * Math.abs(rightDef[0] - that.cx) + /* Math.abs(rightDef[1] - that.cy) + */ 2 * Math.abs(rightDef[2] - that.cz);

            let cpr = rightDistance - leftDistance;
            
            if (leftDistance === rightDistance) {
                cpr = Math.abs(rightDef[1] - that.cy) - Math.abs(leftDef[1] - that.cy);
            }

            return cpr;
        });

        this.removeTasks = this.removeTasks.sort(function (left, right) {
            let leftDef = left[1];
            let rightDef = right[1];

            let leftDistance = 2 * Math.abs(leftDef[0] - that.cx) + Math.abs(leftDef[1] - that.cy) + 2 * Math.abs(leftDef[2] - that.cz);
            let rightDistance = 2 * Math.abs(rightDef[0] - that.cx) + Math.abs(rightDef[1] - that.cy) + 2 * Math.abs(rightDef[2] - that.cz);

            return leftDistance - rightDistance;
        });
    }

    fetchYetAnotherTask() {

        this.dominatingChunkActionPriority = 1 - this.dominatingChunkActionPriority;

        if (this.dominatingChunkActionPriority === 1) {
            let task = this.addTasks.pop();
            if (!task)
                task = this.removeTasks.pop();
            return task;
        }
        else {
            let task = this.removeTasks.pop();
            if (!task)
                task = this.addTasks.pop();
            return task;
        }
    }

    processAsyncChain() {

        if (this.asyncRunning === false)
            return;

        let yetAnotherTaskData = this.fetchYetAnotherTask();

        if (yetAnotherTaskData === undefined) {
            this.asyncRunning = false;
            console.log("[TerrainGridManipulator] terrain updated");
            return;
        }

        let yetAnotherTask = yetAnotherTaskData[1];

        let cx = yetAnotherTask[0];
        let cy = yetAnotherTask[1];
        let cz = yetAnotherTask[2];
        let taskType = yetAnotherTask[3];

        let didNothing = true;

        let isTaskClose = Math.abs(cx - this.cx) + Math.abs(cy - this.cy) + Math.abs(cz - this.cz) <= 5;

        // add, but unless chunk is not too far
        if (taskType === 1) {
            if (!this.terrain.hasChunk(cx, cy, cz)) {
                if (Math.max(Math.abs(this.cx - cx), Math.abs(this.cy - cy), Math.abs(this.cz - cz)) <= CHUNKS_VISIBLE_RADIUS_FAR) {
                    let tc = this.chunkProvider.obtainChunkFrom(cx, cy, cz);
                    this.terrain.uploadChunk(tc);
                    didNothing = isTaskClose;
                }
            }
        }
        
        // delete unless chunk is not close enough
        if (taskType === 0) {
            if (Math.max(Math.abs(this.cx - cx), Math.abs(this.cy - cy), Math.abs(this.cz - cz)) > CHUNKS_VISIBLE_RADIUS) {
                this.terrain.dropChunk(cx, cy, cz);
                didNothing = false;
            }
        }

        this.cancelTaskAt(cx, cy, cz);

        /* Repeat something if this batch of asynchronous update did nothing */
        if (didNothing) {
            this.processAsyncChain();
        }
        else {
            let that = this;

            setTimeout(() => {
                that.processAsyncChain();
            }, 0);
        }         
    }

    kickstartAsyncChain() {

        let that = this;

        if (!this.asyncRunning) {
            this.asyncRunning = true;
            console.log("[TerrainGridManipulator] terrain update started");
            setTimeout(() => {
                that.processAsyncChain();
            }, 0);
        }
    }

    reloadRequiredTerrainParts(characterPosition) {
        
        let cx = Math.floor(characterPosition[0] / CHUNK_SIZE);
        let cy = Math.floor(characterPosition[1] / CHUNK_SIZE);
        let cz = Math.floor(characterPosition[2] / CHUNK_SIZE);

        //this.processAsyncChain();

        /* TODO: trigger update only if gone way too far - like, 2 chunks */
        if (this.smallCx === cx && this.smallCy === cy && this.smallCz === cz) {
            return;
        }

        let step = Math.max(Math.abs(this.cx - cx), Math.abs(this.cy - cy), Math.abs(this.cz - cz));
        let isBigStep = step > 3 || step === undefined || isNaN(step);
            
        this.smallCx = cx;
        this.smallCy = cy;
        this.smallCz = cz;

        if (isBigStep) {
            this.cx = cx;
            this.cy = cy;
            this.cz = cz;
        }

        /* Update visibility, mark for deletion if too far away */
        if (isBigStep) {
            for (const [key, value] of Object.entries(this.terrain.chunkMap)) {
                let terr = this.terrain;

                let dist = Math.max(Math.abs(value.cx - cx), Math.abs(value.cy - cy), Math.abs(value.cz - cz));
                        
                if (dist > CHUNKS_VISIBLE_RADIUS_FAR) {
                    this.askAsyncDeleteAt(value.cx, value.cy, value.cz);
                }

                terr.updateVisibility(value.cx, value.cy, value.cz, dist <= CHUNKS_VISIBLE_RADIUS_PSEUDO_FAR);
                
            }
        }

        let selectedRadius = (isBigStep ? CHUNKS_VISIBLE_RADIUS : 3);
        let selectedHeightRadius = (isBigStep ? CHUNKS_VISIBLE_HEIGHT_RADIUS : 3);

        /* Upload if too close */
        for (let i = -selectedRadius + cx; i < selectedRadius + cx; ++i) {
            for (let j = -selectedHeightRadius + cy; j < selectedHeightRadius + cy; ++j) {
                for (let k = -selectedRadius + cz; k < selectedRadius + cz; ++k) {
                    this.askAsyncUploadAt(i, j, k);
                }
            }
        }
        
        this.refreshTaskPriority();

        for (let i = 0; i < 3; ++i)
            this.processAsyncChain();   // at this point, also update task priority          
        
    }

};

/* ========================= ENGINE-DEPENDENT CLASSES =============================================================== */

/* TODO: make it disposable */
class TerrainRendererForBabylon {

    constructor(scene, material) {
        this.scene = scene;
        this.material = material;
        this.lookupMeshes = new Object();
    }

    deleteChunk(x, y, z) {
        
        let mesh = this.lookupMeshes[[x, y, z]];

        if (mesh !== null && mesh !== undefined) {
            mesh.dispose();
            delete this.lookupMeshes[[x, y, z]];
        }
    }

    updateChunk(x, y, z, sgb) {

        if (sgb.positions.length === 0) {
            this.deleteChunk(x, y, z);
            return;
        }

        let scene = this.scene;

        let customMesh = this.lookupMeshes[[x, y, z]];

        if (customMesh === null || customMesh === undefined)
            customMesh = new BABYLON.Mesh("custom", scene);
        
        let vertexData = new BABYLON.VertexData();
        
        vertexData.positions = sgb.positions;
        vertexData.indices = sgb.indices;
        vertexData.normals = sgb.normals;
        vertexData.uvs = sgb.uvs;

        vertexData.applyToMesh(customMesh);

        customMesh.material = this.material;

        this.lookupMeshes[[x, y, z]] = customMesh;

    }

    updateVisibility(x, y, z, visible) {
        let mesh = this.lookupMeshes[[x, y, z]];

        if (mesh) {
            mesh.isVisible = visible;
        }
    }

};



/* ========================= UNPROCESSED CODE ===================================================================== */


const canvas = document.getElementById("renderCanvas"); // Get the canvas element
const engine = new BABYLON.Engine(canvas, true); // Generate the BABYLON 3D engine

let myCustomMesh;

let myMaterial;

let myCamera;

// Add your code here matching the playground format
const createScene = function () {

    const scene = new BABYLON.Scene(engine);
    scene.clearColor = new BABYLON.Color3(0.3, 0.3, 0.6);

    //const camera = new BABYLON.ArcRotateCamera("camera", -Math.PI / 2, Math.PI / 2.5, 15, new BABYLON.Vector3(0, 0, 0));
    const camera = new BABYLON.UniversalCamera("UniversalCamera", new BABYLON.Vector3(0, -5, -10), scene);

    camera.attachControl(canvas, true);
    camera.position.y = 50;
    //camera.setTarget(new BABYLON.Vector3(0, 0, 0));

    myCamera = camera;
    
    const light = new BABYLON.HemisphericLight("light", new BABYLON.Vector3(1, 1, 0));

    //-------------------------------------------------------------------------------------

    const groundMat2 = new BABYLON.StandardMaterial("groundMat");
    groundMat2.diffuseColor = new BABYLON.Color3(0, 1, 0);

    var img = new Image(512, 512);
    img.src = 'images/trava.png';
    
    const boxMat2 = new BABYLON.StandardMaterial("boxMat");

    let myTextr = new BABYLON.Texture("images/atlas.png", scene, false, undefined);
    
    //myTextr.updateSamplingMode(BABYLON.Texture.NEAREST_SAMPLINGMODE);
    myTextr.wrapR = BABYLON.CLAMP_ADDRESSMODE;
    
    myTextr.anisotropicFilteringLevel = 1;

    boxMat2.diffuseTexture = myTextr;

    //--------------------------------------------------------------------------------------

    myMaterial = boxMat2;
    
    //var customMesh = new BABYLON.Mesh("custom", scene);

    //myCustomMesh = customMesh;

    //customMesh.material = boxMat2;

    //-------------------------------------------------------------------------------------

    return scene;
};

const scene = createScene(); //Call the createScene function


let terrainRenderBabylon = new TerrainRendererForBabylon(scene, myMaterial);

let terrain = new TerrainGrid(terrainRenderBabylon);

let chunkProvider = new TerrainChunkProvider();

let terrainManip = new TerrainGridManipulator(terrain, chunkProvider);
terrainManip.generate();




// Register a render loop to repeatedly render the scene
engine.runRenderLoop(function () {
    
    terrainManip.reloadRequiredTerrainParts([myCamera.position.x, myCamera.position.y, myCamera.position.z]);

    scene.render();
});

// Watch for browser/canvas resize events
window.addEventListener("resize", function () {
        engine.resize();
});
