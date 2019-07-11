#include "CompressionTool.h"


//zlib
GZStreamAdpater::GZStreamAdpater()
{

}

void GZStreamAdpater::open(std::string& filename)
{
    file.open(filename.c_str());
}

void GZStreamAdpater::close()
{
    file.close();
}

std::ostream &GZStreamAdpater::operator <<(std::string& stream)
{
    return file << stream;
}

std::ostream &GZStreamAdpater::operator<<(std::ostream & stream)
{
    return file << stream.rdbuf();
}


//LZ4
#ifdef LZ4
LZ4Adapter::LZ4Adapter() : os(file)
{

}

void LZ4Adapter::open(std::string & filepath)
{
    file = std::ofstream(filepath);
}

void LZ4Adapter::close()
{
    os.close();
}

std::ostream &LZ4Adapter::operator <<(std::string& stream)
{
    return os << stream;
}

std::ostream &LZ4Adapter::operator<<(std::ostream & stream)
{
    return os << stream.rdbuf();
}
#endif

//Uncompressed
UncompressedAdapter::UncompressedAdapter()
{

}

void UncompressedAdapter::open(std::string & filepath)
{
    file = std::ofstream(filepath, std::ofstream::out);
}

void UncompressedAdapter::close()
{
    file.close();
}

std::ostream &UncompressedAdapter::operator<<(std::string & stream)
{
    return file << stream;
}

std::ostream &UncompressedAdapter::operator<<(std::ostream & stream)
{
    return file << stream.rdbuf();
}


//ZStandard
#ifdef ZSTD
ZSTDAdapter::ZSTDAdapter(int level) : os(file,level)
{
}

void ZSTDAdapter::open(std::string & filepath)
{
    file = std::ofstream(filepath);
}

void ZSTDAdapter::close()
{
}

std::ostream &ZSTDAdapter::operator<<(std::string & stream)
{
    return os << stream;
}

std::ostream &ZSTDAdapter::operator<<(std::ostream & stream)
{
    return os << stream.rdbuf();
}
#endif
