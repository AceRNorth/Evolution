#include <vector>
#include <algorithm>
#include "GDRelease.h"
#include "Patch.h"
#include "globals.h"
#include "random.h"
#include <iostream> // for error messages

/**
 * Releases the gene drive mosquitoes into the simulation area if it is time for release.
 * @param[in] day 			simulation day
 * @param[in, out] sites	vector of all Patch objects
 */
void GDRelease::release_gene_drive(int day, std::vector<Patch*> &sites) 
{
    if (is_release_time(day)) {
        auto rel_sites = select_driver_sites(day, sites);
        put_driver_sites(rel_sites);
    }
    if (day==mutant_time) {

    int num_rel_sites = std::min(int(sites.size()), num_mut_sites);
    std::vector<int> rel_patches; // patches in which to release the mutants
	while (rel_patches.size() < num_rel_sites) {
		int rel_pat = random_discrete(0, sites.size() - 1);
		auto is_unique = (rel_patches.end() == std::find(rel_patches.begin(), rel_patches.end(), rel_pat));
		if (is_unique) { 
			rel_patches.push_back(rel_pat);
		}
	}
	    //int rel_pat = random_discrete(0, sites.size() - 1);
	for(int i=0;i<rel_patches.size();i++)
	{
		int rel_pat=rel_patches[i];
		sites[rel_pat]->add_mutants(num_mutants);
	}
//	std::cerr<<" mutant release in patch "<<rel_pat<<"   "<<num_mutants<<std::endl;
    }
}

/**
 * Determines if the day is a chosen release time.
 * @param[in] day simulation day
 * @return As you would expect.
 */
bool GDRelease::is_release_time(int day)
{
	return (std::find(release_times.begin(), release_times.end(), day) == release_times.end()) ? false : true;
}

/**
 * Adds drive heterozygous (WD) male mosquitoes to the release sites.
 * @details Relevant parameters include the number of male driver mosquitoes to be released at each site ``num_driver_M``.
 * @param[in, out] rel_sites release sites
 */
void GDRelease::put_driver_sites(std::vector<Patch*>& rel_sites) 
{
    for (const auto& rel_pat : rel_sites) 
    	{
		rel_pat->add_driver_M(num_driver_M);
	}
}

/**
 * RandomGDRelease constructor.
 * @param[in] params gene drive release parameters
 */
RandomGDRelease::RandomGDRelease(ReleaseParams* params): GDRelease(params->num_driver_M, params->release_times)
{
	num_driver_sites = params->num_driver_sites;
//	num_mutants = params->num_mutants;
//	mutant_time = params->mutant_time;
}

/**
 * Selects random release sites.
 * @details Relevant parameters include the number of driver sites.
 * @param[in] day	simulation day
 * @param[in] sites vector of all Patch objects
 * @return The chosen release sites.
 * @see RandomGDRelease::num_driver_sites
 */
std::vector<Patch*> RandomGDRelease::select_driver_sites(int day, const std::vector<Patch*> &sites)
{
    int num_rel_sites = std::min(int(sites.size()), num_driver_sites);
    std::vector<Patch*> rel_patches; // patches in which to release the gene drive
	while (rel_patches.size() < num_rel_sites) {
		int rel_pat = random_discrete(0, sites.size() - 1);

		// only pick unique sites within the central area to release the gene drive
		auto is_unique = (rel_patches.end() == std::find(rel_patches.begin(), rel_patches.end(), sites[rel_pat]));
		if (is_unique) { 
			rel_patches.push_back(sites[rel_pat]);
		}
	}
	return rel_patches;
}

/**
 * SchedGDRelease constructor. 
 * @param[in] params 	gene drive release parameters
 * @param[in] rel_sites chosen release sites, indices relative to Model::get_sites()
 * @param[in] sites		vector of all Patch objects
 */
SchedGDRelease::SchedGDRelease(ReleaseParams* params, std::vector<int> rel_sites, std::vector<Patch*> &sites): GDRelease(params->num_driver_M, params->release_times) 
{	
    for (const auto& s : rel_sites) { // convert indices relative to sites vector into patch pointers
        release_sites.push_back(sites.at(s));
    }
}

/**
 * Selects the release sites corresponding to the selected release time.
 * @param[in] day	simulation day
 * @param[in] sites vector of all Patch objects
 * @return The release sites.
 */
std::vector<Patch*> SchedGDRelease::select_driver_sites(int day, const std::vector<Patch*> &sites)
{
	return release_sites;
}
