#ifndef GENERALMETAPOP_GDRELEASE_H
#define GENERALMETAPOP_GDRELEASE_H

#include <vector>
#include "Params.h"
#include "Patch.h"

class Patch;

/**
 * Gene drive release base class. Defines the basic implementation of gene drive release in the model.
 * @details GDRelease classes implement the selection mechanism of the release sites from available patches and the release of gene drive mosquitoes into the Patch objects at the chosen release times. The number of drive heterozygous male mosquitoes released into each release site at each release time is determined by ``num_driver_M``. 
 */
class GDRelease {
public:
	/**
	 * GDRelease constructor.
	 * @param[in] num_driver_M 	number of drive heterozygous (WD) male mosquitoes per release
	 * @param[in] rel_times 	days on which the gene drive mosquitoes will be released
	 */
	GDRelease(int num_driver_M, std::vector<int> rel_times): num_driver_M(num_driver_M), release_times(rel_times) {}
	~GDRelease() {};
	void release_gene_drive(int day, std::vector<Patch*> &sites);
protected:
	int num_driver_M;/**< Number of drive heterozygous (WD) male mosquitoes per release. */
	std::vector<int> release_times; /**< Days on which the gene drive mosquitoes will be released. */ 
	bool is_release_time(int day);
    virtual std::vector<Patch*> select_driver_sites(int day, const std::vector<Patch*> &sites) = 0;
	void put_driver_sites(std::vector<Patch*>& rel_sites);
};

/**
 * Implements gene drive release for randomised release sites. 
 * @details Release sites are randomised upon each release time, with the number selected on each release dependent on ``num_driver_sites``.
 * @see ReleaseParams::num_driver_sites
 */
class RandomGDRelease: public GDRelease {
public:
    RandomGDRelease(ReleaseParams* params);
private:
    int num_driver_sites; /**< Number of gene drive release sites per release. */ 
   // int mutant_time; 
   // int num_mutants; 
    std::vector<Patch*> select_driver_sites(int day, const std::vector<Patch*> &sites) override;
};

/**
 * Implements gene drive release for pre-selected release sites.
 * @details Release sites are kept the same for each release, using the same pre-selected sites from the available Patches. 
 */
class SchedGDRelease: public GDRelease {
public:
	SchedGDRelease(ReleaseParams* params, std::vector<int> rel_sites, std::vector<Patch*> &sites);
private:
	std::vector<Patch*> release_sites; /**< Gene drive release sites. */ 
    std::vector<Patch*> select_driver_sites(int day, const std::vector<Patch*> &sites) override;
};

#endif //GENERALMETAPOP_GDRELEASE_H
