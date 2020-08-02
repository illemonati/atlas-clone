cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(bamtoolsLib CXX C)

add_library(${PROJECT_NAME}
        api/api_global.h
        api/BamAlgorithms.h
        api/BamAlignment.cpp api/BamAlignment.h
        api/BamAux.h
        api/BamConstants.h
        api/BamIndex.h
        api/BamMultiReader.cpp api/BamMultiReader.h
        api/BamReader.cpp api/BamReader.h
        api/BamWriter.cpp api/BamWriter.h
        api/IBamIODevice.h
        api/SamConstants.h
        api/SamHeader.cpp api/SamHeader.h
        api/SamProgramChain.cpp api/SamProgramChain.h
        api/SamProgram.cpp api/SamProgram.h
        api/SamReadGroup.cpp
        api/SamReadGroupDictionary.cpp api/SamReadGroupDictionary.h
        api/SamReadGroup.h
        api/SamSequence.cpp api/SamSequence.h
        api/SamSequenceDictionary.cpp api/SamSequenceDictionary.h

        api/algorithms/Sort.h

        api/internal/bam/BamHeader_p.cpp api/internal/bam/BamHeader_p.h
        api/internal/bam/BamMultiMerger_p.h
        api/internal/bam/BamMultiReader_p.cpp api/internal/bam/BamMultiReader_p.h
        api/internal/bam/BamRandomAccessController_p.cpp api/internal/bam/BamRandomAccessController_p.h
        api/internal/bam/BamReader_p.cpp api/internal/bam/BamReader_p.h
        api/internal/bam/BamWriter_p.cpp api/internal/bam/BamWriter_p.h

        api/internal/index/BamIndexFactory_p.cpp api/internal/index/BamIndexFactory_p.h
        api/internal/index/BamStandardIndex_p.cpp api/internal/index/BamStandardIndex_p.h
        api/internal/index/BamToolsIndex_p.cpp api/internal/index/BamToolsIndex_p.h

        api/internal/io/BamDeviceFactory_p.cpp api/internal/io/BamDeviceFactory_p.h
        api/internal/io/BamFile_p.cpp api/internal/io/BamFile_p.h
        api/internal/io/BamFtp_p.cpp api/internal/io/BamFtp_p.h
        api/internal/io/BamHttp_p.cpp api/internal/io/BamHttp_p.h
        api/internal/io/BamPipe_p.cpp api/internal/io/BamPipe_p.h
        api/internal/io/BgzfStream_p.cpp api/internal/io/BgzfStream_p.h
        api/internal/io/ByteArray_p.cpp api/internal/io/ByteArray_p.h
        api/internal/io/HostAddress_p.cpp api/internal/io/HostAddress_p.h
        api/internal/io/HostInfo_p.cpp api/internal/io/HostInfo_p.h
        api/internal/io/HttpHeader_p.cpp api/internal/io/HttpHeader_p.h
        api/internal/io/ILocalIODevice_p.cpp api/internal/io/ILocalIODevice_p.h
        api/internal/io/NetUnix_p.h
        api/internal/io/NetWin_p.h
        api/internal/io/RollingBuffer_p.cpp api/internal/io/RollingBuffer_p.h
        api/internal/io/TcpSocketEngine_p.cpp api/internal/io/TcpSocketEngine_p.h
        api/internal/io/TcpSocketEngine_unix_p.cpp
        api/internal/io/TcpSocket_p.cpp api/internal/io/TcpSocket_p.h

        api/internal/sam/SamFormatParser_p.cpp
        api/internal/sam/SamFormatParser_p.h
        api/internal/sam/SamFormatPrinter_p.cpp
        api/internal/sam/SamFormatPrinter_p.h
        api/internal/sam/SamHeaderValidator_p.cpp
        api/internal/sam/SamHeaderValidator_p.h
        api/internal/sam/SamHeaderVersion_p.h

        api/internal/utils/BamException_p.cpp api/internal/utils/BamException_p.h

        utils/bamtools_fasta.cpp utils/bamtools_fasta.h
        utils/bamtools_filter_engine.h
        utils/bamtools_filter_properties.h
        utils/bamtools_filter_ruleparser.h
        utils/bamtools_options.cpp utils/bamtools_options.h
        utils/bamtools_pileup_engine.cpp utils/bamtools_pileup_engine.h
        utils/bamtools_utilities.cpp utils/bamtools_utilities.h
        utils/bamtools_variant.h
        utils/utils_global.h
        )

include_directories(".")