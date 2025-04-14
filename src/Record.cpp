#include <filesystem>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "Record.h"
#include "constants.h"

/**
 * Record constructor.
 * Creates LocalData, Totals and CoordinateList output .txt files.
 * @details Creates a subdirectory for output files in the current directory.
 * @param[in] rec_params 	recording parameters
 * @param[in] rep 			initial repetition label for the given set of runs
 */
Record::Record(RecordParams* rec_params, int rep) 
{
	rec_start = rec_params->rec_start;
	rec_end = rec_params->rec_end;
	rec_interval_global = rec_params->rec_interval_global;
	rec_interval_local = rec_params->rec_interval_local;
	rec_sites_freq = rec_params->rec_sites_freq;
	set_label = rec_params->set_label;
	rep_label = rep;

	// create folder for destination of output files 
	if (!std::filesystem::exists("./output_files")) {
		std::filesystem::create_directory("output_files");
	}
	std::filesystem::current_path("./output_files");
	
	os1 << "LocalData" << set_label << "run" << rep_label << ".txt"; 
	local_data.open(os1.str());
	os2 << "Totals" << set_label << "run" << rep_label << ".txt";
	global_data.open(os2.str());
	os3 << "CoordinateList" << set_label << "run" << rep_label << ".txt";
	coord_list.open(os3.str());

	local_data << "Male populations of each allele at each site\n";

	if(constants::num_gen==6)
	{
	local_data << "Day" << "\t" << "Site" << "\t" << "W" << "\t" << "D" << "\t"  << "R" << std::endl;
	}
	if(constants::num_gen==10)
	{
	local_data << "Day" << "\t" << "Site" << "\t" << "W" << "\t" << "D" << "\t"  << "R1" << "\t" << "R2" << std::endl;
	}

	if(constants::num_gen==18)
	{
	local_data << "Day" << "\t" << "Site" << "\t" << "W" << "\t" << "D" << "\t"  << "R" << "\t" << "A" << "\t" << "B" << std::endl;
	}
	//local_data << "Male populations of each genotype at each site\n";
	//local_data << "Day" << "\t" << "Site" << "\t" << "WWAA" << "\t" << "WDAA" << "\t" << "DDAA" << "\t" << "WRAA" << "\t" << "RRAA" << "\t" << "DRAA" << "\t" << "WWAB" << "\t" << "WDAB" << "\t" << "DDAB" << "\t" << "WRAB" << "\t" << "RRAB" << "\t" << "DRAB" << "\t" << "WWBB" << "\t" << "WDBB" << "\t" << "DDBB" << "\t" << "WRBB" << "\t" << "RRBB" << "\t" << "DRBB" << std::endl;

	global_data << "Total males of each genotype\n";
	if(constants::num_gen==6)
	{
	global_data << "Day" << "\t" << "WW" << "\t" << "WD" << "\t" << "DD" << "\t" << "WR" << "\t" << "RR" << "\t" << "DR" << std::endl;
	}

	if(constants::num_gen==10)
	{
	global_data << "Day" << "\t" << "WW" << "\t" << "WD" << "\t" << "WR1" << "\t" << "WR2" << "\t" << "DD" << "\t" << "DR1" << "\t" << "DR2" << "\t" << "R1R1" << "\t" << "R1R2" << "\t" << "R2R2" << "\t" << std::endl;
	}

	if(constants::num_gen==18)
	{
	global_data << "Day" << "\t" << "WWAA" << "\t" << "WDAA" << "\t" << "DDAA" << "\t" << "WRAA" << "\t" << "RRAA" << "\t" << "DRAA" << "\t" << "WWAB" << "\t" << "WDAB" << "\t" << "DDAB" << "\t" << "WRAB" << "\t" << "RRAB" << "\t" << "DRAB" << "\t" << "WWBB" << "\t" << "WDBB" << "\t" << "DDBB" << "\t" << "WRBB" << "\t" << "RRBB" << "\t" << "DRBB" << std::endl;
	}
	coord_list << "Coordinate list of the sites\n";
	coord_list << "Site" << "\t" << "x" << "\t" << "y" << std::endl;
}

/**
 * Record destructor.
 * Resets the current filepath so the output_files directory can be found in the next set of runs. 
 */
Record::~Record()
{
	std::filesystem::current_path("..");
}

/**
 * @brief Records the coordinates of the population sites. 
 * @details Relevant parameters include the fraction of sites to collect data for.
 * @param[in] sites vector of all Patch objects
 * @see InputParams::rec_sites_freq
 */
void Record::record_coords(const std::vector<Patch*> &sites) 
{
	const auto default_precision{std::cout.precision()};
	constexpr auto max_precision{std::numeric_limits<double>::digits10 + 1};
	for (int pat=0; pat < sites.size(); pat += rec_sites_freq) {
		auto coords = sites[pat]->get_coords();
		coord_list << pat+1 << "\t" << std::setprecision(max_precision) << coords.x << "\t" << coords.y << std::endl;
	}
	coord_list << std::setprecision(default_precision);
}

/**
 * @brief Records the total numbers of male mosquitoes for the given day, divided by genotype. 
 * @details The totals are assumed to be across all sites. 
 * @param[in] day 		simulation day
 * @param[in] tot_M_gen total number of males divided by genotype
 * @see Model::calculate_tot_M_gen(), Patch::get_M()
 */
void Record::record_global(int day, const std::array<long long int, constants::num_gen> &tot_M_gen)
{
	global_data << day;
	for (const auto& m_gen : tot_M_gen) {
		global_data << "\t" << m_gen;
	}
	global_data << std::endl;
}

/**
 * @brief Outputs the total numbers of juvenile (J), male (M), virgin female (V) and mated female (F) mosquitoes for the given day.
 * @details The totals are assumed to be across all sites, and over all genotypes and age groups. 
 * @param[in] day 	simulation day
 * @param[in] tot_J	total number of juveniles
 * @param[in] tot_M	total number of males
 * @param[in] tot_V	total number of virgin (unmated) females
 * @param[in] tot_F	total number of mated females
 * @see Patch
 */
void Record::output_totals(int day, long long int tot_J, long long int tot_M, long long int tot_V, long long int tot_F)
{
	if (day == 0) {
		std::cout << "\n" << "rep " << rep_label << "\n";
		std::cout << "day" << "   " << "total J" << "   " << "total M" << "   " << "total V" << "   " << "total F" << "\n";
	}
	std::cout << day << "     " << tot_J << "   " << tot_M << "   " << tot_V << "   " << tot_F << std::endl;
}

/**
 * @brief Records the number of males at each site for the given day.
 * @details The number of males at each site is divided by genotype. Relevant parameters include the fraction of sites to collect data for.
 * @param[in] day 	simulation day
 * @param[in] sites vector of all Patch objects
 * @see InputParams::rec_sites_freq
 */
void Record::record_local(int day, const std::vector<Patch*> &sites) 
{
	int W,D,R,A,B,R1,R2;
	for (int pat=0; pat < sites.size(); pat += rec_sites_freq) {
		local_data << day << "\t" << pat+1;

	const auto& m_gen_vector = sites[pat]->get_M();
	if(m_gen_vector.size()==6)
	{
	W= 2*m_gen_vector[0]+m_gen_vector[1]+m_gen_vector[3];
	D= 2*m_gen_vector[2]+m_gen_vector[1]+m_gen_vector[5];
	R= 2*m_gen_vector[4]+m_gen_vector[3]+m_gen_vector[5];
	local_data << "\t" << W << "\t" << D << "\t" << R;
	}

	if(m_gen_vector.size()==10)
	{
	W= 2*m_gen_vector[0]+m_gen_vector[1]+ +m_gen_vector[2] +m_gen_vector[3];
	D= 2*m_gen_vector[4]+m_gen_vector[1]+ +m_gen_vector[5] +m_gen_vector[6];
	R1= 2*m_gen_vector[7]+m_gen_vector[2]+ +m_gen_vector[5] +m_gen_vector[8];
	R2= 2*m_gen_vector[9]+m_gen_vector[3]+ +m_gen_vector[6] +m_gen_vector[8];
	local_data << "\t" << W << "\t" << D << "\t" << R1 << "\t" << R2;
	}

	if(m_gen_vector.size()==18)
	{
	W= 2*m_gen_vector[0]+m_gen_vector[1]+m_gen_vector[3]+ 2*m_gen_vector[6]+m_gen_vector[7]+m_gen_vector[9]+ 2*m_gen_vector[12]+m_gen_vector[13]+m_gen_vector[15];
	D= 2*m_gen_vector[2]+m_gen_vector[1]+m_gen_vector[5]+ 2*m_gen_vector[8]+m_gen_vector[7]+m_gen_vector[11]+ 2*m_gen_vector[14]+m_gen_vector[13]+m_gen_vector[17];
	R= 2*m_gen_vector[4]+m_gen_vector[3]+m_gen_vector[5]+ 2*m_gen_vector[10]+m_gen_vector[9]+m_gen_vector[11]+ 2*m_gen_vector[16]+m_gen_vector[15]+m_gen_vector[17];


	A= 2*(m_gen_vector[0]+m_gen_vector[1]+m_gen_vector[2]+ m_gen_vector[3]+m_gen_vector[4]+m_gen_vector[5])+ (m_gen_vector[6]+m_gen_vector[7]+m_gen_vector[8]+m_gen_vector[9]+m_gen_vector[10]+m_gen_vector[11]);
	B= 2*(m_gen_vector[12]+m_gen_vector[13]+m_gen_vector[14]+ m_gen_vector[15]+m_gen_vector[16]+m_gen_vector[17])+ (m_gen_vector[6]+m_gen_vector[7]+m_gen_vector[8]+m_gen_vector[9]+m_gen_vector[10]+m_gen_vector[11]);
	local_data << "\t" << W << "\t" << D << "\t" << R << "\t" << A << "\t" << B;
	}

	//	for (const auto& m_gen : sites[pat]->get_M()) {
	//		local_data << "\t" << m_gen;
	//	}
	local_data << std::endl;
	}
}

/**
 * @brief Determines if it is time to record global data.
 * @details The number of males at each site is divided by genotype.  
 * @param[in] day 	simulation day
 * @return As you would expect.
 * @see Record::record_global()
 */
bool Record::is_rec_global_time(int day)
{
	return day % rec_interval_global == 0;
}

/**
 * @brief Determines if it is time to record local data.
 * @note The initialisation day (day 0) will always be recorded, and the recording window will be inclusive of the start and end times. 
 * @details Other relevant parameters include the local recording interval. 
 * @param[in] day 	simulation day
 * @return As you would expect.
 * @see Record::record_local(), InputParams::rec_start, InputParams::rec_end, InputParams::rec_interval_local
 */
bool Record::is_rec_local_time(int day) 
{
	return (day == 0) || (day >= rec_start && day <= rec_end && day % rec_interval_local == 0);
}
