# Set Name of project and language
project(atlas LANGUAGES CXX)

# Set cmake version
cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 14)

add_subdirectory(commonutilities)

# include sub-directories (that don't contain a CMakeLists.txt)
include_directories(GLF)
include_directories(PopulationTools)
include_directories(Recalibration)
include_directories(Simulations)
include_directories(Tests)
include_directories(Vcf)
include_directories(bamtools)

# Create executable using the specified src
add_executable(${PROJECT_NAME}
        allTasks.h
        atlas.cpp
        QualityTables.cpp
        QualityTables.h
        TAlignment.cpp
        TAlignment.h
        TAlignmentMerger.cpp
        TAlignmentMerger.h
        TAlignmentParser.cpp
        TAlignmentParser.h
        TAllelicDepthCounts.cpp
        TAllelicDepthCounts.h
        TBase.cpp
        TBase.h
        TBed.cpp
        TBed.h
        TBedReader.cpp
        TBedReader.h
        TCaller.cpp
        TCaller.h
        TChromosomes.cpp
        TChromosomes.h
        TDistributionOfCounts.h
        TFastaBuffer.cpp
        TFastaBuffer.h
        TGenome.cpp
        TGenome.h
        TGenotypeMap.h
        TGenotypePrior.h
        TPostMortemDamage.cpp
        TPostMortemDamage.h
        TQualityMap.h
        TReadGroups.cpp
        TReadGroups.h
        TReadList.h
        TSite.cpp
        TSite.h
        TSiteSubset.cpp
        TSiteSubset.h
        TSoftClipping.cpp
        TSoftClipping.h
        TThetaEstimator.cpp
        TThetaEstimatorData.cpp
        TThetaEstimatorData.h
        TThetaEstimator.h
        TWindow.cpp
        TWindow.h

        PopulationTools/TAlleleCountEstimator.cpp
        PopulationTools/TAlleleCountEstimator.h
        PopulationTools/TAlleleCountFileFormat.cpp
        PopulationTools/TAlleleCountFileFormat.h
        PopulationTools/TAlleleFrequencyEstimator.cpp
        PopulationTools/TAlleleFrequencyEstimator.h
        PopulationTools/TDistanceEstimator.cpp
        PopulationTools/TDistanceEstimator.h
        PopulationTools/TGenotypeFrequencies.cpp
        PopulationTools/TGenotypeFrequencies.h
        PopulationTools/TInbreedingEstimator.cpp
        PopulationTools/TInbreedingEstimator.h
        PopulationTools/TMajorMinor.cpp
        PopulationTools/TMajorMinor.h
        PopulationTools/TPolymorhicWindowIdentifier.cpp
        PopulationTools/TPolymorhicWindowIdentifier.h
        PopulationTools/TPopulationLikelihoodLocus.cpp
        PopulationTools/TPopulationLikelihoodLocus.h
        PopulationTools/TPopulationLikelihoods.cpp
        PopulationTools/TPopulationLikelihoods.h

        Recalibration/TRecalibrationBQSR.cpp
        Recalibration/TRecalibrationBQSR.h
        Recalibration/TRecalibration.cpp
        Recalibration/TRecalibrationEMAuxiliaryTools.cpp
        Recalibration/TRecalibrationEMAuxiliaryTools.h
        Recalibration/TRecalibrationEMEstimator.cpp
        Recalibration/TRecalibrationEMEstimator.h
        Recalibration/TRecalibrationEMModel.cpp
        Recalibration/TRecalibrationEMModel.h
        Recalibration/TRecalibration.h

        GLF/TGLF.cpp GLF/TGLF.h
        GLF/TGlfMultiReader.cpp GLF/TGlfMultiReader.h

        Simulations/SFS.cpp
        Simulations/SFS.h
        Simulations/TSimulatorAuxiliaryTools.cpp
        Simulations/TSimulatorAuxiliaryTools.h
        Simulations/TSimulator.cpp
        Simulations/TSimulator.h
        Simulations/TSimulatorQualityTransformation.cpp
        Simulations/TSimulatorQualityTransformation.h
        Simulations/TSimulatorRead.cpp
        Simulations/TSimulatorRead.h
        Simulations/TSimulatorReadLength.cpp
        Simulations/TSimulatorReadLength.h

        Vcf/TVCFCompare.cpp
        Vcf/TVCFCompare.h
        Vcf/TVcfConverter.cpp
        Vcf/TVcfConverter.h
        Vcf/TVcfDiagnostics.cpp
        Vcf/TVcfDiagnostics.h
        Vcf/TVCFFields.h
        Vcf/TVcfFile.cpp
        Vcf/TVcfFile.h
        Vcf/TVcfParser.cpp
        Vcf/TVcfParser.h

        Tests/TAtlasTest.cpp
        Tests/TAtlasTestFilter.cpp
        Tests/TAtlasTestFilter.h
        Tests/TAtlasTest.h
        Tests/TAtlasTestMergePairs.cpp
        Tests/TAtlasTestMergePairs.h
        Tests/TAtlasTestPMD.cpp
        Tests/TAtlasTestPMD.h
        Tests/TAtlasTestRecalibration.cpp
        Tests/TAtlasTestRecalibration.h
        Tests/TVcfTest.cpp
        Tests/TVcfTest.h

        bamtools/api/api_global.h
        bamtools/api/BamAlgorithms.h
        bamtools/api/BamAlignment.cpp bamtools/api/BamAlignment.h
        bamtools/api/BamAux.h
        bamtools/api/BamConstants.h
        bamtools/api/BamIndex.h
        bamtools/api/BamMultiReader.cpp bamtools/api/BamMultiReader.h
        bamtools/api/BamReader.cpp bamtools/api/BamReader.h
        bamtools/api/BamWriter.cpp bamtools/api/BamWriter.h
        bamtools/api/IBamIODevice.h
        bamtools/api/SamConstants.h
        bamtools/api/SamHeader.cpp bamtools/api/SamHeader.h
        bamtools/api/SamProgramChain.cpp bamtools/api/SamProgramChain.h
        bamtools/api/SamProgram.cpp bamtools/api/SamProgram.h
        bamtools/api/SamReadGroup.cpp
        bamtools/api/SamReadGroupDictionary.cpp bamtools/api/SamReadGroupDictionary.h
        bamtools/api/SamReadGroup.h
        bamtools/api/SamSequence.cpp bamtools/api/SamSequence.h
        bamtools/api/SamSequenceDictionary.cpp bamtools/api/SamSequenceDictionary.h

        bamtools/api/algorithms/Sort.h

        bamtools/api/internal/bam/BamHeader_p.cpp bamtools/api/internal/bam/BamHeader_p.h
        bamtools/api/internal/bam/BamMultiMerger_p.h
        bamtools/api/internal/bam/BamMultiReader_p.cpp bamtools/api/internal/bam/BamMultiReader_p.h
        bamtools/api/internal/bam/BamRandomAccessController_p.cpp bamtools/api/internal/bam/BamRandomAccessController_p.h
        bamtools/api/internal/bam/BamReader_p.cpp bamtools/api/internal/bam/BamReader_p.h
        bamtools/api/internal/bam/BamWriter_p.cpp bamtools/api/internal/bam/BamWriter_p.h

        bamtools/api/internal/index/BamIndexFactory_p.cpp bamtools/api/internal/index/BamIndexFactory_p.h
        bamtools/api/internal/index/BamStandardIndex_p.cpp bamtools/api/internal/index/BamStandardIndex_p.h
        bamtools/api/internal/index/BamToolsIndex_p.cpp bamtools/api/internal/index/BamToolsIndex_p.h

        bamtools/api/internal/io/BamDeviceFactory_p.cpp bamtools/api/internal/io/BamDeviceFactory_p.h
        bamtools/api/internal/io/BamFile_p.cpp bamtools/api/internal/io/BamFile_p.h
        bamtools/api/internal/io/BamFtp_p.cpp bamtools/api/internal/io/BamFtp_p.h
        bamtools/api/internal/io/BamHttp_p.cpp bamtools/api/internal/io/BamHttp_p.h
        bamtools/api/internal/io/BamPipe_p.cpp bamtools/api/internal/io/BamPipe_p.h
        bamtools/api/internal/io/BgzfStream_p.cpp bamtools/api/internal/io/BgzfStream_p.h
        bamtools/api/internal/io/ByteArray_p.cpp bamtools/api/internal/io/ByteArray_p.h
        bamtools/api/internal/io/HostAddress_p.cpp bamtools/api/internal/io/HostAddress_p.h
        bamtools/api/internal/io/HostInfo_p.cpp bamtools/api/internal/io/HostInfo_p.h
        bamtools/api/internal/io/HttpHeader_p.cpp bamtools/api/internal/io/HttpHeader_p.h
        bamtools/api/internal/io/ILocalIODevice_p.cpp bamtools/api/internal/io/ILocalIODevice_p.h
        bamtools/api/internal/io/NetUnix_p.h
        bamtools/api/internal/io/NetWin_p.h
        bamtools/api/internal/io/RollingBuffer_p.cpp bamtools/api/internal/io/RollingBuffer_p.h
        bamtools/api/internal/io/TcpSocketEngine_p.cpp bamtools/api/internal/io/TcpSocketEngine_p.h
        bamtools/api/internal/io/TcpSocketEngine_unix_p.cpp
        bamtools/api/internal/io/TcpSocket_p.cpp bamtools/api/internal/io/TcpSocket_p.h

        bamtools/api/internal/sam/SamFormatParser_p.cpp
        bamtools/api/internal/sam/SamFormatParser_p.h
        bamtools/api/internal/sam/SamFormatPrinter_p.cpp
        bamtools/api/internal/sam/SamFormatPrinter_p.h
        bamtools/api/internal/sam/SamHeaderValidator_p.cpp
        bamtools/api/internal/sam/SamHeaderValidator_p.h
        bamtools/api/internal/sam/SamHeaderVersion_p.h

        bamtools/api/internal/utils/BamException_p.cpp bamtools/api/internal/utils/BamException_p.h

        bamtools/utils/bamtools_fasta.cpp bamtools/utils/bamtools_fasta.h
        bamtools/utils/bamtools_filter_engine.h
        bamtools/utils/bamtools_filter_properties.h
        bamtools/utils/bamtools_filter_ruleparser.h
        bamtools/utils/bamtools_options.cpp bamtools/utils/bamtools_options.h
        bamtools/utils/bamtools_pileup_engine.cpp bamtools/utils/bamtools_pileup_engine.h
        bamtools/utils/bamtools_utilities.cpp bamtools/utils/bamtools_utilities.h
        bamtools/utils/bamtools_variant.h
        bamtools/utils/utils_global.h)

# load packages
find_package(Armadillo REQUIRED)
find_package(ZLIB REQUIRED)
target_include_directories(${PROJECT_NAME} PUBLIC ${ARMADILLO_INCLUDE_DIRS}
        "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/commonutilities")

# Define required c++ standard to C++14
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)

# link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC commonutilities ${ARMADILLO_LIBRARIES} ${ZLIB_LIBRARIES})