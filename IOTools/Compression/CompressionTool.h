#ifndef COMPRESSIONTOOL_H
#define COMPRESSIONTOOL_H

#include <string>
#include <cstdio>
#include <ostream>
#include <sstream>
#include "gzstream.h"

#ifdef LZ4
#include "lz4_stream.h"
#endif

#ifdef ZSTD
#include "zstd_stream.h"
#endif

class CompressionTool
{
public:
    virtual void open(std::string&)=0;
    virtual void close()=0;
    virtual std::ostream& operator<<(std::string&)=0;
    virtual std::ostream& operator<<(std::ostream&)=0;

    virtual ~CompressionTool(){}
};

class GZStreamAdpater : public CompressionTool{
private:
    gz::ogzstream file;
    // CompressionTool interface
public:
    GZStreamAdpater();

    void open(std::string&);
    void close();
    std::ostream &operator<<(std::string&);
    std::ostream &operator<<(std::ostream&);

    ~GZStreamAdpater(){}
};

#ifdef LZ4
class LZ4Adapter : public CompressionTool{

    std::ofstream file;
    lz4_stream::ostream os;

    // CompressionTool interface
public:
    LZ4Adapter();

    void open(std::string&);
    void close();
    std::ostream &operator<<(std::string&);
    std::ostream &operator<<(std::ostream&);

    ~LZ4Adapter(){}
};
#endif


#ifdef ZSTD
class ZSTDAdapter : public CompressionTool{

    std::ofstream file;
    compressed_streams::ZstdOStream os;

    // CompressionTool interface
public:
    ZSTDAdapter(int);

    void open(std::string&);
    void close();
    std::ostream &operator<<(std::string&);
    std::ostream &operator<<(std::ostream&);

    ~ZSTDAdapter(){}
};
#endif

class UncompressedAdapter : public CompressionTool{

    std::ofstream file;

    // CompressionTool interface
public:
    UncompressedAdapter();

    void open(std::string&);
    void close();
    std::ostream &operator<<(std::string&);
    std::ostream &operator<<(std::ostream&);

    ~UncompressedAdapter(){}
};

#endif // COMPRESSIONTOOL_H
