#ifndef TALLELES_H_
#define TALLELES_H_

#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/GenomePositions/TGenomeWindow.h"
#include "genometools/GenotypeTypes.h"
#include "coretools/Containers/TView.h"

namespace SiteSubset {

struct TAllele {
	genometools::TGenomePosition position;
	genometools::Base ref;
	genometools::Base alt;

	TAllele(size_t RefID, size_t Pos, genometools::Base Ref, genometools::Base Alt)
		: position(RefID, Pos), ref(Ref), alt(Alt) {}

	friend bool operator<(const TAllele& lhs, const TAllele& rhs) noexcept {return lhs.position < rhs.position;}
	friend bool operator<(const TAllele& lhs, const genometools::TGenomePosition& rhs) noexcept {return lhs.position < rhs;}
	friend bool operator<(const genometools::TGenomePosition& lhs, const TAllele& rhs) noexcept {return lhs < rhs.position;}
};

enum class Morphic {Mono, Poly};

class TAlleles {
	using Sites = coretools::TConstView<TAllele>;
	std::vector<std::string> _chromosomes;
	std::vector<TAllele> _allSites;
	std::vector<Sites> _sites;

	size_t _nChrSites = 0;
		
public:
	using iterator       = Sites::iterator;
	using const_iterator = Sites::const_iterator;

	TAlleles() = default;
	TAlleles(std::string_view Filename, Morphic M = Morphic::Poly) { parse(Filename, M); }
	TAlleles(std::string_view Filename, const genometools::TChromosomes &Chromosomes, Morphic M = Morphic::Poly) {parse(Filename, Chromosomes, M);}
	TAlleles(std::string_view Filename, coretools::TConstView<std::string> Chromosomes, Morphic M = Morphic::Poly) {parse(Filename, Chromosomes, M);}

	void parse(std::string_view Filename, Morphic M = Morphic::Poly); 
	void parse(std::string_view Filename, const genometools::TChromosomes &Chromosomes, Morphic M = Morphic::Poly);
	void parse(std::string_view Filename, coretools::TConstView<std::string> Chromosomes, Morphic M = Morphic::Poly);


	bool empty() const noexcept { return _allSites.empty();}
	operator bool() const noexcept {return !empty(); }

	size_t NChrSites() const noexcept { return _nChrSites; }
	size_t NChr() const noexcept { return _chromosomes.size(); }

	size_t refID(std::string_view Name) const noexcept {
		return std::distance(_chromosomes.cbegin(), std::find(_chromosomes.cbegin(), _chromosomes.cend(), Name));
	}
	size_t size() const noexcept { return _allSites.size(); }
	const std::string &name(size_t RefID) const { return _chromosomes[RefID]; }


	const Sites &sites(size_t RefID) const noexcept { return _sites[RefID]; }
	const Sites &sites(std::string_view Name) const noexcept { return sites(refID(Name)); }
	template<typename Key> const Sites &operator[](const Key& K) const noexcept { return sites(K); }

	const_iterator begin() const noexcept { return _allSites.data(); }
	const_iterator end() const noexcept { return _allSites.data() + _allSites.size(); }

	template<typename Key> size_t size(const Key& K) const noexcept { return window(K).size(); }
	template<typename Key> bool empty(const Key& K) const noexcept { return window(K).empty(); }

	const_iterator begin(size_t K) const noexcept { return sites(K).begin(); }
	const_iterator begin(std::string_view K) const noexcept { return sites(K).begin(); }
	const_iterator begin(const genometools::TGenomeWindow &Window) const noexcept {
		const auto &ws = sites(Window.refID());

		auto begin = std::lower_bound(ws.begin(), ws.end(), Window.from());
		if (begin != ws.begin()) begin--;

		for (auto it = begin; it != ws.end(); ++it) {
			if (Window.within(it->position)) return it;
			if (it->position > Window.to()) break;
		}
		return end();
	}

	bool overlaps(const genometools::TGenomeWindow &Window) const noexcept { return (begin(Window) != end()); }
	bool overlaps(const genometools::TGenomePosition &Position) const noexcept { return (begin(genometools::TGenomeWindow{Position, 1}) != end()); }

};
}

#endif
