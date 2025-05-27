#ifndef TBGZWRITER_H_
#define TBGZWRITER_H_

#include "coretools/Main/TError.h"
#include "htslib/bgzf.h"

#include "coretools/Files/TWriter.h"

namespace GLF {

class TBGzWriter final : public coretools::TWriter {
	BGZF *_file;

	void _write(const void *buffer, size_t size, size_t count) override {
		if (bgzf_write(_file, buffer, size * count) == 0) { throw coretools::TDevError("Was not able to write to gz file!"); };
	};
	int64_t _tell() const override { return bgzf_tell(_file); };

public:
	TBGzWriter(std::string_view Filename, const char *Mode = "w")
		: TWriter(Filename), _file(bgzf_open(name().c_str(), Mode)) {
		if (!_file) { UERROR("Was not able to create file ", name(), ". Does the path exist?"); }
	}
	~TBGzWriter() { bgzf_close(_file); }

	TBGzWriter(const TBGzWriter &)            = delete;
	TBGzWriter &operator=(const TBGzWriter &) = delete;

	TBGzWriter(TBGzWriter &&)                 = default;
	TBGzWriter &operator=(TBGzWriter &&)      = default;
};

} // namespace coretools

#endif
