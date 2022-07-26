//Tamar Schlick lab All rights reserved.

#include "trajectory.h"
#include "utilities.h"

long int trajectory::number_of_frames = 0;
long int trajectory::number_of_analyzed_frames = 0;

using namespace std;

int trajectory::fill_core_particles(string file_name){

    ifstream core_file(file_name.c_str());
    int i = 0;

    if (!core_file) {
        printf("The core file %s has not been opened\n", file_name.c_str());
        return -1;
    }
    else {
        while ((core_file >> core_particles.coordinates[i][0] >> core_particles.coordinates[i][1] >> core_particles.coordinates[i][2]) && (i<300)){
            core_particles.coordinates[i][0] = core_particles.coordinates[i][0] / 10;
            core_particles.coordinates[i][1] = core_particles.coordinates[i][1] / 10;
            core_particles.coordinates[i][2] = core_particles.coordinates[i][2] / 10;
            i++;
        }
        core_file.close();
        cout << endl << "The core file has been successfully read!" << endl;
        return 0;
    }
}



int trajectory::fill_LH_data(string file_name){
    string line, segment;
    ifstream LH_file(file_name.c_str());
    stringstream ss;
    int i = 0, p = 0;

    if (!LH_file) {
        printf("The LH file %s has not been opened\n", file_name.c_str());
        return -1;
    }
    else{

        while (getline(LH_file, line)){

            if (line[0] != '#') {

                ss << line;

                i = 0;
                while (std::getline(ss, segment, ' ')){
                    if (segment.length() > 0){
                        i++;
                        if (i == 5){

                            this->LH_ev_hh.push_back(atof(segment.c_str()));

                            if (p<6)

                                this->GH_ev_hh.push_back(atof(segment.c_str()));
                            else

                                this->CTD_ev_hh.push_back(atof(segment.c_str()));

                            p++;
                        }

                        if (i == 6){

                            this->LH_ev_hc.push_back(atof(segment.c_str()));

                            if (p<6)

                                this->GH_ev_hc.push_back(atof(segment.c_str()));
                            else

                                this->CTD_ev_hc.push_back(atof(segment.c_str()));
                        }

                        if (i == 7){

                            this->LH_ev_hl.push_back(atof(segment.c_str()));

                            if (p<6)

                                this->GH_ev_hl.push_back(atof(segment.c_str()));
                            else

                                this->CTD_ev_hl.push_back(atof(segment.c_str()));
                        }


                        if (i == 8){

                            this->LH_ev_ht.push_back(atof(segment.c_str()));

                            if (p<6)

                                this->GH_ev_ht.push_back(atof(segment.c_str()));
                            else

                                this->CTD_ev_ht.push_back(atof(segment.c_str()));
                            ss.str("");
                            break;
                        }

                    }
                }
            }
        }

        LH_file.close();
        cout << "The LH file has been successfully read!" << endl << endl;
        return 0;
    }

}


double trajectory::calculate_radius_of_gyration(vector < row > LHc){

    double mean_x, mean_y, mean_z;
    vector< double > x_diff;
    vector< double > y_diff;
    vector< double > z_diff;
    double radius_of_gyration;
    int i;

    mean_x = 0; mean_y = 0; mean_z = 0;
    for (i = 0; i < this->LH_beads; i++){
        mean_x = mean_x + LHc[i].x;
        mean_y = mean_y + LHc[i].y;
        mean_z = mean_z + LHc[i].z;
    }

    mean_x = mean_x / LHc.size();
    mean_y = mean_y / LHc.size();
    mean_z = mean_z / LHc.size();

    for (i = 0; i < this->LH_beads; i++){
        x_diff.push_back(LHc[i].x - mean_x);
        y_diff.push_back(LHc[i].y - mean_y);
        z_diff.push_back(LHc[i].z - mean_z);
    }

    radius_of_gyration = 0;
    for (i = 0; i < this->LH_beads; i++)
        radius_of_gyration = radius_of_gyration + x_diff[i] * x_diff[i] + y_diff[i] * y_diff[i] + z_diff[i] * z_diff[i];

    radius_of_gyration = sqrt(radius_of_gyration / LHc.size());
    vector<double>().swap(x_diff); vector<double>().swap(y_diff); vector<double>().swap(z_diff);

    return radius_of_gyration;
}

int trajectory::parse_trajectory(string file_name){

    //    double radius_of_gyration;
    int line_number = 0;
    long int frame_number = 0;
    char text[256];
    char* pEnd;
    int counter_DNA = 0, counter_n = 0;
    row coords;
    row_orientation o_coords;


    cout << "path/filename = " << file_name << "\n";
    printf("Setup :\n");
    printf("cores : %d\n", this->cores);
    printf("DNA beads : %d\n", this->DNA_beads);
    printf("Tail beads : %d\n", this->tail_beads);
    printf("LH beads : %d\n", this->LH_beads);
    printf("lines_per_frame = %d\n", this->lines_per_frame);

    FILE* simulation_file = fopen(file_name.c_str(), "r");

    if (simulation_file){

        cout << "Parsing file!" << endl;

        while (fgets(text, sizeof(text), simulation_file))
        {

            line_number++;


            if ((line_number <= this->lines_per_frame) & ((frame_number + 1) >= this->starting_frame)) {

                //cout << endl << "Filling data!" << endl;

                // Read nucleosome positions and orientations
                if ((line_number >= this->starting_line_for_nucleosome) & (line_number <= this->ending_line_for_nucleosome)){

                    counter_n++;
                    if (counter_n == 1){
                        coords.x = strtod(text, &pEnd);
                        coords.y = strtod(pEnd, &pEnd);
                        coords.z = strtod(pEnd, NULL);
                        this->nucleosome_c.push_back(coords);
                    }
                    else{
                        o_coords.a = strtod(text, &pEnd);
                        o_coords.b = strtod(pEnd, &pEnd);
                        o_coords.c = strtod(pEnd, NULL);
                        this->nucleosome_orientations.push_back(o_coords);
                    }

                    if (counter_n == 4)
                        counter_n = 0;
                   
                }

              
                // Read DNA positions and orientations
                if ((line_number >= this->starting_line_for_DNA) & (line_number <= this->ending_line_for_DNA)){
                    coords.x = strtod(text, &pEnd);
                    coords.y = strtod(pEnd, &pEnd);
                    coords.z = strtod(pEnd, NULL);
                    counter_DNA++;

                    if (counter_DNA == 1){
                        coords.x = strtod(text, &pEnd);
                        coords.y = strtod(pEnd, &pEnd);
                        coords.z = strtod(pEnd, NULL);
                        this->DNA_c.push_back(coords);
                    }
                    else{
                        o_coords.a = strtod(text, &pEnd);
                        o_coords.b = strtod(pEnd, &pEnd);
                        o_coords.c = strtod(pEnd, NULL);
                        this->DNA_orientations.push_back(o_coords);
                    }

                    if (counter_DNA == 4)
                        counter_DNA = 0;
                }


                // Read tails positions
                if ((line_number >= this->starting_line_for_tails) & (line_number <= this->ending_line_for_tails)){
                    coords.x = strtod(text, &pEnd);
                    coords.y = strtod(pEnd, &pEnd);
                    coords.z = strtod(pEnd, NULL);
                    this->tail_c.push_back(coords);
                }


                // Read LH positions
                if ((line_number >= this->starting_line_for_LH) & (line_number <= this->ending_line_for_LH)) {
                    coords.x = strtod(text, &pEnd);
                    coords.y = strtod(pEnd, &pEnd);
                    coords.z = strtod(pEnd, NULL);
                    this->LH_c.push_back(coords);

                    if ((line_number - this->starting_line_for_LH + 1) % this->LH_beads == 0){ // calculates rg after each LH                    
                            } // if ((line_number - this->starting_line_for_LH + 1) % 28 == 0)
                    //} // if ((frame_number + 1) >= this->starting_frame){/
                } //if ((line_number >= this->starting_line_for_LH) & (line_number <= this->ending_line_for_LH)) {

            }  // if ((line_number <= this->lines_per_frame) & ((frame_number + 1) >= this->starting_frame)) {





            if (line_number == this->lines_per_frame){
                line_number = 0;
                frame_number++;
                this->number_of_frames++;



                if ((frame_number) >= this->starting_frame){

                    number_of_analyzed_frames++;

                    this->totdata++;

                  
                    this->calculate_core_core_interactions();
                    
                    vector< row >().swap(this->LH_c);             // erase LH_c
                    vector< row >().swap(this->nucleosome_c);     // erase nucleosome_c
                    vector< row_orientation >().swap(this->nucleosome_orientations);
                    vector< row >().swap(this->DNA_c);            // erase DNA_c
                    vector< row_orientation >().swap(this->DNA_orientations); // erase DNA_orientations
                    vector< row >().swap(this->tail_c);           // erase tail_c
                }

                
                if ((frame_number + 1) < this->starting_frame)
                    cout << ".";
                else
                    cout << "o";
            }
        }

        fclose(simulation_file);
    }

    
    cout << endl << endl;

    return 0;
}


int trajectory::set_parameters(int cores, int DNA_beads, int LH_beads, int tail_beads, int starting_frame, bool save_data){

    row test;

    this->cores = cores;
    this->DNA_beads = DNA_beads;
    this->LH_beads = LH_beads;
    this->tail_beads = tail_beads;
    this->lines_per_frame = cores * 4 + DNA_beads * 4 * cores + cores*tail_beads + cores*LH_beads;


    this->starting_line_for_nucleosome = 1;
    this->ending_line_for_nucleosome = cores * 4;

    this->starting_line_for_DNA = cores * 4 + 1;
    this->ending_line_for_DNA = cores * 4 + cores * DNA_beads * 4;

    this->starting_line_for_tails = cores * 4 + cores * DNA_beads * 4 + 1;
    this->ending_line_for_tails = cores * 4 + cores * DNA_beads * 4 + cores*tail_beads;
    

    this->starting_line_for_LH = cores * 4 + DNA_beads * 4 * cores + cores*tail_beads + 1;
    this->ending_line_for_LH = cores * 4 + DNA_beads * 4 * cores + cores*tail_beads + cores*LH_beads;
   
    this->starting_frame = starting_frame;
    this->save_data = save_data;


       interCC.resize(cores, std::vector<int>(cores, 0));

    this->totdata = 0;

    cout << "Matrix set!" << endl;

    cout << "nucleosome_c.size() = ";
    cout << nucleosome_c.size() << endl;

    test.x = 0.0;
    test.y = 0.1;
    test.z = 0.2;

    nucleosome_c.push_back(test);

    return 0;
}



int trajectory::save_to_file(){
    if (this->save_data){
        if (this->rg.size() > 0){
            char file_name[256];
            sprintf(file_name, "rg_cores_%d_beads_%d.txt", this->cores, this->DNA_beads);

            ofstream outputfile(file_name);
            if (outputfile){
                outputfile.precision(5);
                for (vector< double >::const_iterator lhrg = this->rg.begin(); lhrg != this->rg.end(); ++lhrg){
                    outputfile << *lhrg << endl;
                }
                outputfile.close();
                return 0;
            }
            else
                return -1;
        }
        return 0;
    }

    return 0;
}


int trajectory::calculate_LH_LH_interactions(){
    int i, j, k, l;
    int core_counter_i = 0, core_counter_k = 0;
    double radius_i = 0, radius_k = 0;
    double dr, dx, dy, dz;
    bool skip_LH_bead = false;
    int count_LH_LH_contacts;

    vector < int > LH_contacts(this->LH_beads, 0);

    for (i = 0; i < this->cores - 1; i++){  // analyzed linker (belonging to core i)

        count_LH_LH_contacts = 0;
        vector< int >().swap(LH_contacts);
        LH_contacts.resize(this->LH_beads, 0);

        for (j = 0; j < this->LH_beads; j++){ // its beads

            skip_LH_bead = false;

            core_counter_i = i * this->LH_beads + j;

            radius_i = this->LH_ev_hh[j];

            for (k = i + 1; k < this->cores; k++){
                if (k != i){  // left if I want to change something (do not want to do reformating)
                    for (l = 0; l < this->LH_beads; l++){
                        core_counter_k = k*this->LH_beads + l;
                        radius_k = this->LH_ev_hh[l];
                        dx = LH_c[core_counter_i].x - LH_c[core_counter_k].x;
                        dy = LH_c[core_counter_i].y - LH_c[core_counter_k].y;
                        dz = LH_c[core_counter_i].z - LH_c[core_counter_k].z;

                        dr = sqrt(dx*dx + dy*dy + dz*dz);

                        if (dr <= (radius_i + radius_k)){
                            LH_contacts[j] = 1;
                            skip_LH_bead = true;
                            break;
                        }

                    } // for (l = 0; l < this->LH_beads; l++){
                } // if (k != i){          
                if (skip_LH_bead)
                    break;
            } // for (k = i+1; k < this->cores; k++){    

        }  //  for (j = 0; j < this->LH_beads; j++){ // its beads


        this->frequency_of_LH_LH_interactions.push_back(arithmetic_mean(LH_contacts));

    }  // for (i = 0; i < this->cores-1; i++){  // analyzed linker (belonging to core i)

    //vector< int >().swap(LH_contacts);

    return 0;
}




int trajectory::calculate_LH_DNA_interactions(){
    int i, j, k, l;
    int core_counter = 0, DNA_counter = 0;
    double radius_i = 0, radius_k = 0;
    double dr, dx, dy, dz;
    bool skip_LH = false;
    bool interacts_with_parental;
    bool interacts_with_nonparental;

    vector < int > LH_parental_DNA_contacts(this->LH_beads, 0);
    vector < int > LH_nonparental_DNA_contacts(this->LH_beads, 0);


    for (i = 0; i < this->cores; i++){  // analyzed linker histone belonging to core i


        LH_parental_DNA_contacts.resize(this->LH_beads, 0);
        LH_nonparental_DNA_contacts.resize(this->LH_beads, 0);

        for (j = 0; j < this->LH_beads; j++){ // its beads

            interacts_with_parental = false;
            interacts_with_nonparental = false;

            core_counter = i * this->LH_beads + j;
            radius_i = this->LH_ev_hl[j];


            for (k = 0; k < (this->cores + 1); k++){
                for (l = (k)*this->DNA_beads; l < (k + 1)*this->DNA_beads; l++){
                    DNA_counter = l;
                    radius_k = DNA_ev;

                    dx = LH_c[core_counter].x - DNA_c[DNA_counter].x;
                    dy = LH_c[core_counter].y - DNA_c[DNA_counter].y;
                    dz = LH_c[core_counter].z - DNA_c[DNA_counter].z;

                    dr = sqrt(dx*dx + dy*dy + dz*dz);

                    if (dr <= (radius_i + radius_k)){
                        if ((i == k) || ((i + 1) == k)){
                            if (!interacts_with_parental){
                                LH_parental_DNA_contacts[j] = 1;
                                interacts_with_parental = true;
                            }

                        }
                        else{
                            LH_nonparental_DNA_contacts[j] = 1;
                            interacts_with_nonparental = true;
                        } // if ((i == k) || (i == k + 1)){
                    } // if (dr <= (radius_i + radius_k)){

                    if (interacts_with_nonparental & interacts_with_nonparental)
                        break;

                }  // for (l = (k)*this->DNA_beads; l < (k)*this->DNA_beads + this->DNA_beads; l++){

                if (interacts_with_nonparental & interacts_with_nonparental)
                    break;

            } // for (k = 0; k < (this->cores+1); k++){            

        }  // for (j = 0; j < this->LH_beads; j++){ // its beads

        this->frequency_of_LH_parental_DNA_interactions.push_back(arithmetic_mean(LH_parental_DNA_contacts));
        this->frequency_of_LH_nonparental_DNA_interactions.push_back(arithmetic_mean(LH_nonparental_DNA_contacts));


    }  //  for (i = 0; i < this->cores; i++){  // analyzed linker histone belonging to core i


   

    return 0;
}


int trajectory::calculate_core_coordinates(int core_number){
    int i, j;

  

    for (i = 0; i < 300; i++){

        j = core_number * 3;

      
        local_core_particles.coordinates[i][0] = nucleosome_orientations[j].a*core_particles.coordinates[i][0] +
            nucleosome_orientations[j + 1].a*core_particles.coordinates[i][1] +
            nucleosome_orientations[j + 2].a*core_particles.coordinates[i][2];

        local_core_particles.coordinates[i][1] = nucleosome_orientations[j].b*core_particles.coordinates[i][0] +
            nucleosome_orientations[j + 1].b*core_particles.coordinates[i][1] +
            nucleosome_orientations[j + 2].b*core_particles.coordinates[i][2];

        local_core_particles.coordinates[i][2] = nucleosome_orientations[j].c*core_particles.coordinates[i][0] +
            nucleosome_orientations[j + 1].c*core_particles.coordinates[i][1] +
            nucleosome_orientations[j + 2].c*core_particles.coordinates[i][2];


     

    }

   

    return 0;
}

int trajectory::calculate_LH_core_interactions(){
    int i, j, k, l;
    int LH_counter = 0, DNA_counter = 0;
    double dr, dx, dy, dz;
    double core_distance_x, core_distance_y, core_distance_z;
    double local_LH_bead_x, local_LH_bead_y, local_LH_bead_z;
    double radius_i, radius_k;
    bool interacts_with_parental;
    bool interacts_with_nonparental;

    vector< core > core_vector;
    vector < int > LH_parental_core_contacts(this->LH_beads, 0);
    vector < int > LH_nonparental_core_contacts(this->LH_beads, 0);


    for (i = 0; i < this->cores; i++){
        this->calculate_core_coordinates(i);
        core_vector.push_back(local_core_particles);
    }


    for (i = 0; i < this->cores; i++){  // analyzed linker histone   

        LH_parental_core_contacts.resize(this->LH_beads, 0);
        LH_nonparental_core_contacts.resize(this->LH_beads, 0);

        for (j = 0; j < this->LH_beads; j++){ // the beads of analyzed linker histone

            interacts_with_parental = false;
            interacts_with_nonparental = false;

            radius_i = this->LH_ev_hc[j];


            LH_counter = i * this->LH_beads + j;

            local_LH_bead_x = this->LH_c[LH_counter].x - this->nucleosome_c[i].x;
            local_LH_bead_y = this->LH_c[LH_counter].y - this->nucleosome_c[i].y;
            local_LH_bead_z = this->LH_c[LH_counter].z - this->nucleosome_c[i].z;


            for (k = 0; k < this->cores; k++){ // analyzed core

                core_distance_x = this->nucleosome_c[i].x - this->nucleosome_c[k].x;
                core_distance_y = this->nucleosome_c[i].y - this->nucleosome_c[k].y;
                core_distance_z = this->nucleosome_c[i].z - this->nucleosome_c[k].z;



                for (l = 0; l<300; l++){

                    radius_k = core_ev;

                    dx = local_LH_bead_x - core_vector[k].coordinates[l][0] + core_distance_x;
                    dy = local_LH_bead_y - core_vector[k].coordinates[l][1] + core_distance_y;
                    dz = local_LH_bead_z - core_vector[k].coordinates[l][2] + core_distance_z;

                    dr = sqrt(dx*dx + dy*dy + dz*dz);

                    if (dr <= (radius_i + radius_k)){

                        if (i == k){
                            LH_parental_core_contacts[j] = 1;
                            interacts_with_parental = true;
                        }
                        else{
                            LH_nonparental_core_contacts[j] = 1;
                            interacts_with_nonparental = true;
                        }
                    }

                    if (interacts_with_parental & interacts_with_nonparental)
                        break;
                }  // for (l = 0; l>300; l++){

                if (interacts_with_parental & interacts_with_nonparental)
                    break;
            } // for (k = 0; k < this->cores; k++){ // analyzed core


        }

        this->frequency_of_LH_parental_core_interactions.push_back(arithmetic_mean(LH_parental_core_contacts));
        this->frequency_of_LH_nonparental_core_interactions.push_back(arithmetic_mean(LH_nonparental_core_contacts));

    }

    

    return 0;
}






int trajectory::calculate_core_core_interactions(){
    int i, j, k, l;
    double dr, dx, dy, dz;
    double cc_distance;
    vector< core > core_vector;
    bool exit_loop = false;

    for (i = 0; i < this->cores; i++){
        this->calculate_core_coordinates(i);
        core_vector.push_back(local_core_particles);
    }

    for (i = 0; i < this->cores-1; i++){  // analyzes tails belonging to the core i
        cout << endl << "Core(i) = " << i << ", Core(k) = ";
        
        for (k = i + 1; k < this->cores; k++){  // analyzes core k

            exit_loop = false;
            
            cc_distance = sqrt((this->nucleosome_c[k].x - this->nucleosome_c[i].x)*(this->nucleosome_c[k].x - this->nucleosome_c[i].x) +
                (this->nucleosome_c[k].y - this->nucleosome_c[i].y)*(this->nucleosome_c[k].y - this->nucleosome_c[i].y) +
                (this->nucleosome_c[k].z - this->nucleosome_c[i].z)*(this->nucleosome_c[k].z - this->nucleosome_c[i].z));

            if(cc_distance <= CORE_CORE_CUTOFF){
                if (i != k){   

                    cout << k << " ";

                    for (j = (i)*this->tail_beads; j < (i + 1)*this->tail_beads; j++){
                        for (l = 0; l < NUMBER_OF_CORE_PARTICLES; l++){                    
                            dx = this->tail_c[j].x - (core_vector[k].coordinates[l][0] + this->nucleosome_c[k].x);
                            dy = this->tail_c[j].y - (core_vector[k].coordinates[l][1] + this->nucleosome_c[k].y);
                            dz = this->tail_c[j].z - (core_vector[k].coordinates[l][2] + this->nucleosome_c[k].z);
                                
                            
                            dr = dx*dx + dy*dy + dz*dz;
                                
                            if (dr <= cut_off_squared){
                                this->interCC[i][k]++;
                                this->interCC[k][i]++;
                                exit_loop = true;
                                break;
                            }
                        }
                        if (exit_loop)
                            break;
                    }

                 
                    if(!exit_loop){
                        for (j = (k)*this->tail_beads; j < (k + 1)*this->tail_beads; j++){
                            for (l = 0; l < NUMBER_OF_CORE_PARTICLES; l++){
                                dx = this->tail_c[j].x - (core_vector[i].coordinates[l][0] + this->nucleosome_c[i].x);
                                dy = this->tail_c[j].y - (core_vector[i].coordinates[l][1] + this->nucleosome_c[i].y);
                                dz = this->tail_c[j].z - (core_vector[i].coordinates[l][2] + this->nucleosome_c[i].z);

                                dr = dx*dx + dy*dy + dz*dz;

                                if (dr <= cut_off_squared){
                                    this->interCC[k][i]++;
                                    this->interCC[i][k]++;
                                    exit_loop = true;
                                    break;
                                }
                            }
                            if (exit_loop)
                                break;
                        }
                    }

                   

                    if(!exit_loop){
                        if(cc_distance<=15){
                
                            for (j = 0; j < NUMBER_OF_CORE_PARTICLES; j++){
                                for (l = 0; l < NUMBER_OF_CORE_PARTICLES; l++){
                                    dx = (core_vector[i].coordinates[j][0] + this->nucleosome_c[i].x) -
                                        (core_vector[k].coordinates[l][0] + this->nucleosome_c[k].x);

                                    dy = (core_vector[i].coordinates[j][1] + this->nucleosome_c[i].y) -
                                        (core_vector[k].coordinates[l][1] + this->nucleosome_c[k].y);

                                    dz = (core_vector[i].coordinates[j][2] + this->nucleosome_c[i].z) -
                                        (core_vector[k].coordinates[l][2] + this->nucleosome_c[k].z);

                                   dr = dx*dx + dy*dy + dz*dz;

                                    if (dr <= cut_off_squared){
                                        this->interCC[k][i]++;
                                        this->interCC[i][k]++;
                                        exit_loop = true;
                                        break;
                                    }
                                }
                                if (exit_loop)
                                    break;
                            }
                        }
                    }
            
        
                
             } // if (i != k){  


            } // if(cc_distance <= CORE_CORE_CUTOFF){


        }
    }

        return 0;    
}

int trajectory::save_core_core_interactions(){
    int i, j;
    char file_name[100];

    sprintf(file_name, "CC_interactions_cores_%d_beads_%d_f%d.txt", this->cores, this->DNA_beads, this->starting_frame);
    ofstream core_core_interactions_file(file_name);

        
    for (i = 0; i<this->cores; i++)
    {
        for (j = 0; j<this->cores; j++){            
            core_core_interactions_file << setw(10) << this->interCC[i][j];
        }
        core_core_interactions_file << endl;
    }
    core_core_interactions_file.close();

    return 0;
}


int trajectory::calculate_LH_tail_interactions(){
    int i, j, k, l;
    int core_counter = 0, tail_counter = 0;
    double radius_i = 0, radius_k = 0;
    double dr, dx, dy, dz;
    bool skip_LH = false;


    vector < int > LH_parental_H2A1_contacts(this->LH_beads, 0);
    vector < int > LH_nonparental_H2A1_contacts(this->LH_beads, 0);

    vector < int > LH_parental_H2A2_contacts(this->LH_beads, 0);
    vector < int > LH_nonparental_H2A2_contacts(this->LH_beads, 0);

    vector < int > LH_parental_H2B_contacts(this->LH_beads, 0);
    vector < int > LH_nonparental_H2B_contacts(this->LH_beads, 0);

    vector < int > LH_parental_H3_contacts(this->LH_beads, 0);
    vector < int > LH_nonparental_H3_contacts(this->LH_beads, 0);

    vector < int > LH_parental_H4_contacts(this->LH_beads, 0);
    vector < int > LH_nonparental_H4_contacts(this->LH_beads, 0);

    radius_k = CDT_tail;

    for (i = 0; i < this->cores; i++){  // analyzed linker histone belonging to core i


        LH_parental_H2A1_contacts.resize(this->LH_beads, 0);
        LH_nonparental_H2A1_contacts.resize(this->LH_beads, 0);

        LH_parental_H2A2_contacts.resize(this->LH_beads, 0);
        LH_nonparental_H2A2_contacts.resize(this->LH_beads, 0);

        LH_parental_H2B_contacts.resize(this->LH_beads, 0);
        LH_nonparental_H2B_contacts.resize(this->LH_beads, 0);

        LH_parental_H3_contacts.resize(this->LH_beads, 0);
        LH_nonparental_H3_contacts.resize(this->LH_beads, 0);

        LH_parental_H4_contacts.resize(this->LH_beads, 0);
        LH_nonparental_H4_contacts.resize(this->LH_beads, 0);

        for (j = 0; j < this->LH_beads; j++){ // its beads

            radius_i = this->LH_ev_ht[j];
            core_counter = i * this->LH_beads + j;

            for (k = 0; k < this->cores; k++){
                for (l = (k)*this->tail_beads; l < (k + 1)*this->tail_beads; l++){

                    tail_counter = l;

                    dx = LH_c[core_counter].x - tail_c[tail_counter].x;
                    dy = LH_c[core_counter].y - tail_c[tail_counter].y;
                    dz = LH_c[core_counter].z - tail_c[tail_counter].z;
                    dr = sqrt(dx*dx + dy*dy + dz*dz);

                    if (dr < (radius_i + radius_k)){


                        if (((l % this->tail_beads) > 25) & ((l % this->tail_beads) < 34)){ // H2A1
                           

                            if (i == k){
                                LH_parental_H2A1_contacts[j] = 1;
                            }
                            else{
                                LH_nonparental_H2A1_contacts[j] = 1;
                            }
                        } // if (((l % this->tail_beads) < 34) & ((l % this->tail_beads) > 25)){


                        if (((l % this->tail_beads) > 43) & ((l % this->tail_beads) < 50)){ // H2A2
                           

                            if (i == k){
                                LH_parental_H2A2_contacts[j] = 1;
                            }
                            else{
                                LH_nonparental_H2A2_contacts[j] = 1;
                            }
                        } // if (((l % this->tail_beads) < 50) & ((l % this->tail_beads) > 43)){ // H2A2


                        if (((l % this->tail_beads) > 33) & ((l % this->tail_beads) < 44)){ // H2B
                           

                            if (i == k){
                                LH_parental_H2B_contacts[j] = 1;
                            }
                            else{
                                LH_nonparental_H2B_contacts[j] = 1;
                            }
                        } // if (((l % this->tail_beads) > 33) & ((l % this->tail_beads) < 44)){ // H2B



                        if ((l % this->tail_beads) < 16){ // H3
                           
                            if (i == k){
                                LH_parental_H3_contacts[j] = 1;
                            }
                            else{
                                LH_nonparental_H3_contacts[j] = 1;
                            }
                        } // if ((l % this->tail_beads) < 16){ // H3


                        if (((l % this->tail_beads) > 15) & ((l % this->tail_beads) < 26)){ // H4
                            

                            if (i == k){
                                LH_parental_H4_contacts[j] = 1;
                            }
                            else{
                                LH_nonparental_H4_contacts[j] = 1;
                            }
                        } // if (((l % this->tail_beads) > 33) & ((l % this->tail_beads) < 44)){ // H4


                    } //if (dr < (radius_i + radius_k)){

                } // for (l = (k)*this->tail_beads; l < (k + 1)*this->tail_beads; l++){


            } // for (k = 0; k < this->cores; k++){
        } // for (j = 0; j < this->LH_beads; j++){ // its beads

        this->frequency_of_LH_parental_H2A1_interactions.push_back(arithmetic_mean(LH_parental_H2A1_contacts));
        this->frequency_of_LH_nonparental_H2A1_interactions.push_back(arithmetic_mean(LH_nonparental_H2A1_contacts));

        this->frequency_of_LH_parental_H2A2_interactions.push_back(arithmetic_mean(LH_parental_H2A2_contacts));
        this->frequency_of_LH_nonparental_H2A2_interactions.push_back(arithmetic_mean(LH_nonparental_H2A2_contacts));

        this->frequency_of_LH_parental_H2B_interactions.push_back(arithmetic_mean(LH_parental_H2B_contacts));
        this->frequency_of_LH_nonparental_H2B_interactions.push_back(arithmetic_mean(LH_nonparental_H2B_contacts));


        this->frequency_of_LH_parental_H3_interactions.push_back(arithmetic_mean(LH_parental_H3_contacts));
        this->frequency_of_LH_nonparental_H3_interactions.push_back(arithmetic_mean(LH_nonparental_H3_contacts));

        this->frequency_of_LH_parental_H4_interactions.push_back(arithmetic_mean(LH_parental_H4_contacts));
        this->frequency_of_LH_nonparental_H4_interactions.push_back(arithmetic_mean(LH_nonparental_H4_contacts));


    } // for (i = 0; i < this->cores; i++){  // analyzed linker histone belonging to core i

    return 0;
}




int trajectory::calculate_statistics(){

    cout << fixed << showpoint;
    cout.precision(2);
#ifdef CALCULATE_RADIUS_OF_GYRATION
    double average = accumulate(this->rg.begin(), this->rg.end(), 0.0);
    average = average / this->rg.size();

    cout << endl << "-----------------------------------------------------------" << endl;
    cout << fixed << showpoint;
    cout.precision(2);
    cout << "Mean value of the radius of gyration is         " << average << endl;

    double stdev = 0.0;
    for (vector< double >::const_iterator lhrg = this->rg.begin(); lhrg != this->rg.end(); ++lhrg){
        stdev += (*lhrg - average) * (*lhrg - average);
    }

    stdev = sqrt(stdev / (this->rg.size() - 1));
    cout << "Standard deviation of the radius of gyration is " << stdev << endl;
#endif

    cout << "The average frequency of LH-LH contacts is " << arithmetic_mean(this->frequency_of_LH_LH_interactions) << endl;

    cout << endl << endl;
    cout << "The average frequency of parental LH-DNA contacts     is " << arithmetic_mean(this->frequency_of_LH_parental_DNA_interactions) << endl;
    cout << "The average frequency of non-parental LH-DNA contacts is " << arithmetic_mean(this->frequency_of_LH_nonparental_DNA_interactions) << endl;

    cout << endl << endl;
    cout << "average frequency of parental LH-core contacts     is " << arithmetic_mean(this->frequency_of_LH_parental_core_interactions) << endl;
    cout << "average frequency of non-parental LH-core contacts is " << arithmetic_mean(this->frequency_of_LH_nonparental_core_interactions) << endl;
    cout << "----------------------------------------------------------------------------------------" << endl;
    cout << endl << endl;
    cout << "average frequency of parental LH-H2A1 contacts     is " << arithmetic_mean(this->frequency_of_LH_parental_H2A1_interactions) << endl;
    cout << "average frequency of non-parental LH-H2A1 contacts is " << arithmetic_mean(this->frequency_of_LH_nonparental_H2A1_interactions) << endl;

    cout << endl << endl;
    cout << "average frequency of parental LH-H2A2 contacts     is " << arithmetic_mean(this->frequency_of_LH_parental_H2A2_interactions) << endl;
    cout << "average frequency of non-parental LH-H2A2 contacts is " << arithmetic_mean(this->frequency_of_LH_nonparental_H2A2_interactions) << endl;

    cout << endl << endl;
    cout << "average frequency of parental LH-H2B contacts     is " << arithmetic_mean(this->frequency_of_LH_parental_H2B_interactions) << endl;
    cout << "average frequency of non-parental LH-H2B contacts is " << arithmetic_mean(this->frequency_of_LH_nonparental_H2B_interactions) << endl;

    cout << endl << endl;
    cout << "average frequency of parental LH-H3 contacts     is " << arithmetic_mean(this->frequency_of_LH_parental_H3_interactions) << endl;
    cout << "average frequency of non-parental LH-H3 contacts is " << arithmetic_mean(this->frequency_of_LH_nonparental_H3_interactions) << endl;

    cout << endl << endl;
    cout << "average frequency of parental LH-H4 contacts     is " << arithmetic_mean(this->frequency_of_LH_parental_H4_interactions) << endl;
    cout << "average frequency of non-parental LH-H4 contacts is " << arithmetic_mean(this->frequency_of_LH_nonparental_H4_interactions) << endl;


    return 0;
}

