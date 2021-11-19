#include "utilities.hpp"

using namespace mfem;


Array<int> evalBlockOffsets (Array<int>& blockSizes)
{
    Array<int> blockOffsets(blockSizes.Size()+1);

    blockOffsets[0] = 0;
    for (int i=0; i<blockSizes.Size(); i++) {
        blockOffsets[i+1] = blockSizes[i];
    }
    blockOffsets.PartialSum();

    return blockOffsets;
}

Array<int> evalTemporalBlockOffsets(int minLevel, int maxLevel)
{
    auto blockSizes = evalTemporalBlockSizes(minLevel, maxLevel);
    Array<int> blockOffsets(blockSizes.Size()+1);

    blockOffsets[0] = 0;
    for (int i=0; i<blockSizes.Size(); i++) {
        blockOffsets[i+1] = blockSizes[i];
    }
    blockOffsets.PartialSum();

    return blockOffsets;
}

Array<int> evalTemporalBlockSizes(int minLevel, int maxLevel)
{
    return evalHierarchicalBlockSizes1D(minLevel, maxLevel);
}

Array<int> evalHierarchicalBlockSizes1D (int minLevel, int maxLevel)
{
    int numLevels = (maxLevel - minLevel)+1;
    Array<int> blockSizes(numLevels);
    
    // for the minimum level, the hierarchical basis functions
    // are the same as the standard basis functions
    blockSizes[0]
            = 1+static_cast<int>(std::pow(2, minLevel));
    for (int i=1; i<numLevels; i++) {
        blockSizes[i]
                = static_cast<int>(std::pow(2, minLevel+i-1));
    }
    
    return blockSizes;
}

Array<int> evalSpaceTimeBlockSizes(Array<int>& spatialBlockSizes,
                                   Array<int>& temporalBlockSizes)
{
    int numLevels = spatialBlockSizes.Size();
    Array<int> spaceTimeBlockSizes(numLevels);

    for (int i=0; i<numLevels; i++) {
        spaceTimeBlockSizes[i] = temporalBlockSizes[i]
                * spatialBlockSizes[numLevels - i - 1];
    }

    return spaceTimeBlockSizes;
}
