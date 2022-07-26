//Tamar Schlick lab All rights reserved.

#include "trajectory.h"
#include "utilities.h"
//#define DO_NOT_CALCULATE_EVERYTHING
#define CALCULATE_EVERYTHING

long int trajectory::number_of_frames = 0;
long int trajectory::number_of_analyzed_frames = 0;

using namespace std;

const double trajectory::X_1 = 4.5651;
const double trajectory::Y_1 = 1.4833;
const double trajectory::Z_1 = 1.800;

const double trajectory::X_end = 0.0;
const double trajectory::Y_end = -4.800;
const double trajectory::Z_end = -1.800;

double dot_product(row v1, row v2){
    double product = 0;
    product = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;

    return product;

}


row cross_product(row v1, row v2){
    row product;

    product.x = v1.y*v2.z - v1.z*v2.y;
    product.y = -1 * (v1.x*v2.z - v1.z*v2.x);
    product.z = v1.x*v2.y - v1.y*v2.x;

    return product;
}

double normr(row v){

    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);

}

int trajectory::fill_core_particles(string file_name){

    ifstream core_file(file_name.c_str());
    int i = 0;

    if (!core_file) {
        printf("The core file %s has not been opened\n", file_name.c_str());
        return -1;
    }
    else {
        while ((core_file >> core_particles.coordinates[i][0] >> core_particles.coordinates[i][1] >> core_particles.coordinates[i][2]) && (i < 300)){
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

                            if (p < 6)

                                this->GH_ev_hh.push_back(atof(segment.c_str()));
                            else

                                this->CTD_ev_hh.push_back(atof(segment.c_str()));

                            p++;
                        }

                        if (i == 6){

                            this->LH_ev_hc.push_back(atof(segment.c_str()));

                            if (p < 6)

                                this->GH_ev_hc.push_back(atof(segment.c_str()));
                            else

                                this->CTD_ev_hc.push_back(atof(segment.c_str()));
                        }

                        if (i == 7){

                            this->LH_ev_hl.push_back(atof(segment.c_str()));

                            if (p < 6)

                                this->GH_ev_hl.push_back(atof(segment.c_str()));
                            else

                                this->CTD_ev_hl.push_back(atof(segment.c_str()));
                        }


                        if (i == 8){

                            this->LH_ev_ht.push_back(atof(segment.c_str()));

                            if (p < 6)

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


int trajectory::calculate_gyration_tensor(){

    double Sxx, Sxy, Sxz, Syy, Syz, Szz;
    double mean_x, mean_y, mean_z;
    int i, j, LH_start;
    ostringstream file_name;
    ofstream outfile;

    file_name << "gyration_tensor_c" << this->cores << "_b" << this->DNA_beads << ".txt";
    //cout << file_name.str() << endl;

    outfile.open(file_name.str().c_str(), std::ios_base::app);
    outfile << fixed << showpoint;
    outfile << setw(8);
    outfile.precision(3);
    
    

    for (i = 0; i < this->cores; i++){

        if (this->LH_presence[i]>0){

            LH_start = i*this->LH_beads;
            Sxx = 0.0; Sxy = 0.0; Sxz = 0.0; Syy = 0.0; Syz = 0.0; Szz = 0.0;

            mean_x = 0; mean_y = 0; mean_z = 0;
            for (j = LH_start + 6; j < (LH_start + this->LH_beads); j++){
                mean_x = mean_x + LH_c[j].x;
                mean_y = mean_y + LH_c[j].y;
                mean_z = mean_z + LH_c[j].z;
            }

            mean_x = mean_x / (this->LH_beads - 6);
            mean_y = mean_y / (this->LH_beads - 6);
            mean_z = mean_z / (this->LH_beads - 6);

            for (j = LH_start + 6; j < (LH_start + this->LH_beads); j++){
                Sxx = Sxx + (LH_c[j].x - mean_x)*(LH_c[j].x - mean_x);
                Sxy = Sxy + (LH_c[j].x - mean_x)*(LH_c[j].y - mean_y);
                Sxz = Sxz + (LH_c[j].x - mean_x)*(LH_c[j].z - mean_z);

                Syy = Syy + (LH_c[j].y - mean_y)*(LH_c[j].y - mean_y);
                Syz = Syz + (LH_c[j].y - mean_y)*(LH_c[j].z - mean_z);

                Szz = Szz + (LH_c[j].z - mean_z)*(LH_c[j].z - mean_z);
            }
            Sxx = Sxx / (this->LH_beads - 6);
            Sxy = Sxy / (this->LH_beads - 6);
            Sxz = Sxz / (this->LH_beads - 6);
            Syy = Syy / (this->LH_beads - 6);
            Syz = Syz / (this->LH_beads - 6);
            Szz = Szz / (this->LH_beads - 6);

            outfile << setw(12) << Sxx << setw(12) << Sxy << setw(12) << Sxz << endl;
            outfile << setw(12) << Sxy << setw(12) << Syy << setw(12) << Syz << endl;
            outfile << setw(12) << Sxz << setw(12) << Syz << setw(12) << Szz << endl;
        }
    }

    outfile.close();

    return 0;
}


int trajectory::calculate_radius_of_gyration(){

    double mean_x, mean_y, mean_z;
    vector< double > x_diff;
    vector< double > y_diff;
    vector< double > z_diff;
    double radius_of_gyration;
    //double pom_1, pom_2, pom_3;
    int i, j, LH_start;

    for (i = 0; i < this->cores; i++){

        if (this->LH_presence[i]>0){

            LH_start = i*this->LH_beads;

            mean_x = 0; mean_y = 0; mean_z = 0;
            //for (j = LH_start + 6; j < (LH_start + this->LH_beads); j++){
            for (j = LH_start + 6; j < (LH_start + this->simulated_LH_beads); j++){
                mean_x = mean_x + LH_c[j].x;
                mean_y = mean_y + LH_c[j].y;
                mean_z = mean_z + LH_c[j].z;
            }

            mean_x = mean_x / (this->simulated_LH_beads - 6);
            mean_y = mean_y / (this->simulated_LH_beads - 6);
            mean_z = mean_z / (this->simulated_LH_beads - 6);

            //for (j = LH_start + 6; j < (LH_start + this->LH_beads); j++){
            for (j = LH_start + 6; j < (LH_start + this->simulated_LH_beads); j++){
                x_diff.push_back(LH_c[j].x - mean_x);
                y_diff.push_back(LH_c[j].y - mean_y);
                z_diff.push_back(LH_c[j].z - mean_z);
            }

            radius_of_gyration = 0;
            for (j = 0; j < x_diff.size(); j++){
            

                radius_of_gyration = radius_of_gyration + x_diff[j] * x_diff[j] + y_diff[j] * y_diff[j] + z_diff[j] * z_diff[j];
                
            }

            radius_of_gyration = sqrt(radius_of_gyration / (this->LH_beads - 6));

            this->rg.push_back(radius_of_gyration);
            vector<double>().swap(x_diff); vector<double>().swap(y_diff); vector<double>().swap(z_diff);
        }
        
    }

    return 0;
}


int trajectory::calculate_distance_9_23(){
    double distance;
    int i, LH_start;
    for (i = 0; i < this->cores; i++){
        if (this->LH_presence[i]>0){
            LH_start = i*this->LH_beads;
            distance = sqrt((this->LH_c[LH_start + 9].x - this->LH_c[LH_start + 23].x)*
                (this->LH_c[LH_start + 9].x - this->LH_c[LH_start + 23].x) +

                (this->LH_c[LH_start + 9].y - this->LH_c[LH_start + 23].y)*
                (this->LH_c[LH_start + 9].y - this->LH_c[LH_start + 23].y) +

                (this->LH_c[LH_start + 9].z - this->LH_c[LH_start + 23].z)*
                (this->LH_c[LH_start + 9].z - this->LH_c[LH_start + 23].z));

            this->distance_9_23.push_back(distance);
        }
    }
    return 0;
}



int trajectory::calculate_distance_9_27(){
    double distance;
    int i, LH_start;
    int LH_end;

    LH_end = this->simulated_LH_beads - 1;

    for (i = 0; i < this->cores; i++){
        if (this->LH_presence[i]>0){
            LH_start = i*this->LH_beads;
            distance = sqrt((this->LH_c[LH_start + 9].x - this->LH_c[LH_start + LH_end].x)*
                (this->LH_c[LH_start + 9].x - this->LH_c[LH_start + LH_end].x) +

                (this->LH_c[LH_start + 9].y - this->LH_c[LH_start + LH_end].y)*
                (this->LH_c[LH_start + 9].y - this->LH_c[LH_start + LH_end].y) +

                (this->LH_c[LH_start + 9].z - this->LH_c[LH_start + LH_end].z)*
                (this->LH_c[LH_start + 9].z - this->LH_c[LH_start + LH_end].z));

       

            this->distance_9_27.push_back(distance);
        }
    }
    return 0;
}


int trajectory::calculate_distance_23_27(){
    double distance;
    int i, LH_start;
    int LH_end;

    LH_end = this->simulated_LH_beads - 1;


    for (i = 0; i < this->cores; i++){
        if (this->LH_presence[i]>0){
            LH_start = i*this->LH_beads;

            distance = sqrt((this->LH_c[LH_start + 23].x - this->LH_c[LH_start + LH_end].x)*
                (this->LH_c[LH_start + 23].x - this->LH_c[LH_start + LH_end].x) +

                (this->LH_c[LH_start + 23].y - this->LH_c[LH_start + LH_end].y)*
                (this->LH_c[LH_start + 23].y - this->LH_c[LH_start + LH_end].y) +

                (this->LH_c[LH_start + 23].z - this->LH_c[LH_start + LH_end].z)*
                (this->LH_c[LH_start + 23].z - this->LH_c[LH_start + LH_end].z));

     

            this->distance_23_27.push_back(distance);
        }
    }
    return 0;
}



int trajectory::calculate_distance_6_27(){
    double distance;
    int i, LH_start;
    int start_n, end_n;
    int LH_end;

    LH_end = this->simulated_LH_beads - 1;


    if (this->cores > 2){
        start_n = 1;
        end_n = this->cores - 1;
    }
    else{
        start_n = 0;
        end_n = this->cores;

    }

    for (i = start_n; i < end_n; i++){
        if (this->LH_presence[i]>0){

            LH_start = i*this->LH_beads;
          
            distance = sqrt((this->LH_c[LH_start + 6].x - this->LH_c[LH_start + LH_end].x)*
                (this->LH_c[LH_start + 6].x - this->LH_c[LH_start + LH_end].x) +

                (this->LH_c[LH_start + 6].y - this->LH_c[LH_start + LH_end].y)*
                (this->LH_c[LH_start + 6].y - this->LH_c[LH_start + LH_end].y) +

                (this->LH_c[LH_start + 6].z - this->LH_c[LH_start + LH_end].z)*
                (this->LH_c[LH_start + 6].z - this->LH_c[LH_start + LH_end].z));

            this->distance_6_27.push_back(distance);
        }
    }
    return 0;
}



int trajectory::calculate_distance_6_16(){
    double distance;
    int i, LH_start;
    int start_n, end_n;

    if (this->cores > 2){
        start_n = 1;
        end_n = this->cores - 1;
    }
    else{
        start_n = 0;
        end_n = this->cores;

    }

    for (i = start_n; i < end_n; i++){

        if (this->LH_presence[i]>0){
            LH_start = i*this->LH_beads;
            distance = sqrt((this->LH_c[LH_start + 6].x - this->LH_c[LH_start + 16].x)*
                (this->LH_c[LH_start + 6].x - this->LH_c[LH_start + 16].x) +

                (this->LH_c[LH_start + 6].y - this->LH_c[LH_start + 16].y)*
                (this->LH_c[LH_start + 6].y - this->LH_c[LH_start + 16].y) +

                (this->LH_c[LH_start + 6].z - this->LH_c[LH_start + 16].z)*
                (this->LH_c[LH_start + 6].z - this->LH_c[LH_start + 16].z));

            this->distance_6_16.push_back(distance);
        }
    }
    return 0;
}


int trajectory::calculate_distance_16_27(){
    double distance;
    int i, LH_start;
    int start_n, end_n;
    int LH_end;

    LH_end = this->simulated_LH_beads - 1;

    if (this->cores > 2){
        start_n = 1;
        end_n = this->cores - 1;
    }
    else{
        start_n = 0;
        end_n = this->cores;

    }

    for (i = start_n; i < end_n; i++){
        if (this->LH_presence[i]>0){
            //LH_start = i*this->LH_beads;
            LH_start = i*this->LH_beads;
            distance = sqrt((this->LH_c[LH_start + 16].x - this->LH_c[LH_start + LH_end].x)*
                (this->LH_c[LH_start + 16].x - this->LH_c[LH_start + LH_end].x) +

                (this->LH_c[LH_start + 16].y - this->LH_c[LH_start + LH_end].y)*
                (this->LH_c[LH_start + 16].y - this->LH_c[LH_start + LH_end].y) +

                (this->LH_c[LH_start + 16].z - this->LH_c[LH_start + LH_end].z)*
                (this->LH_c[LH_start + 16].z - this->LH_c[LH_start + LH_end].z));

            /*distance = sqrt((this->LH_c[LH_start + 16].x - this->LH_c[LH_start + 27].x)*
                (this->LH_c[LH_start + 16].x - this->LH_c[LH_start + 27].x) +

                (this->LH_c[LH_start + 16].y - this->LH_c[LH_start + 27].y)*
                (this->LH_c[LH_start + 16].y - this->LH_c[LH_start + 27].y) +

                (this->LH_c[LH_start + 16].z - this->LH_c[LH_start + 27].z)*
                (this->LH_c[LH_start + 16].z - this->LH_c[LH_start + 27].z));*/

            this->distance_16_27.push_back(distance);
        }
    }
    return 0;
}



int trajectory::calculate_distance(char what){
    double distance = 0;
    int index = what * 11 + 6 - 1;
    int i, LH_start;
    for (i = 0; i < this->cores; i++){

        if (this->LH_presence[i]>0){

            //LH_start = i*this->LH_beads;
            LH_start = i*this->simulated_LH_beads;
            distance = sqrt((this->LH_c[LH_start + 6].x - this->LH_c[LH_start + index].x)*(this->LH_c[LH_start + 6].x - this->LH_c[LH_start + index].x) +
                (this->LH_c[LH_start + 6].y - this->LH_c[LH_start + index].y)*(this->LH_c[LH_start + 6].y - this->LH_c[LH_start + index].y) +
                (this->LH_c[LH_start + 6].z - this->LH_c[LH_start + index].z)*(this->LH_c[LH_start + 6].z - this->LH_c[LH_start + index].z));

            if (what == 1)
                this->distance_1_11.push_back(distance);
            else
                this->distance_1_22.push_back(distance);
        }
    }

    return 0;
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
    //int i, j;


    cout << "path/filename = " << file_name << "\n";
    printf("Setup :\n");
    printf("cores : %d\n", this->cores);
    printf("DNA beads : %d\n", this->DNA_beads);
    printf("Tail beads : %d\n", this->tail_beads);
    printf("LH beads : %d\n", this->LH_beads);
    printf("lines_per_frame = %d\n", this->lines_per_frame);

    FILE* simulation_file = fopen(file_name.c_str(), "r");

    if (simulation_file){

        while (fgets(text, sizeof(text), simulation_file))
        {

            line_number++;

            if ((line_number <= this->lines_per_frame) & ((frame_number + 1) >= this->starting_frame)) {

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
                        //radius_of_gyration = calculate_radius_of_gyration(LH_c);
                        //rg.push_back(radius_of_gyration);
                        //vector< row >().swap(LH_c);  // erase LH_c
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

                   
                    /*this->calculate_LH_core_interactions();
                    this->calculate_DNA_bending_angles();
                    this->calculate_angle_between_cores();*/

                   
                    this->calculate_radius_of_gyration();
                    


#ifdef CALCULATE_EVERYTHING

                    this->calculate_LH_DNA_interactions();
                    this->calculate_LH_core_interactions();

                    this->calculate_1_11_22_angle();
                    this->calculate_angle_between_DNAs();
                    this->calculate_angle_between_cores();

                    this->calculate_DNA_bending_angles();
                    this->calculate_LH_LH_interactions();
                    this->calculate_LH_DNA_interactions();
                    this->calculate_LH_core_interactions();
                    this->calculate_LH_tail_interactions();
                                  
                    
                    this->calculate_radius_of_gyration();
                   // this->calculate_core_core_distance();

                    

                    this->distance_between_neighboring_LHs(27, 27);
                    this->distance_between_neighboring_LHs(24, 27);
                    this->distance_between_neighboring_LHs(21, 27);
                    this->distance_between_neighboring_LHs(16, 27);
                    this->distance_between_neighboring_LHs(11, 27);
                    this->distance_between_neighboring_LHs(5, 27);
                                        
                    this->distance_between_neighboring_LHs(16, 21);
                    this->distance_between_neighboring_LHs(11, 21);
                    this->distance_between_neighboring_LHs(5, 21);

                    this->distance_between_neighboring_LHs(5, 16);

                    this->distance_between_neighboring_LHs(26, 26);
                    this->distance_between_neighboring_LHs(25, 25);
                    this->distance_between_neighboring_LHs(24, 24);


                    this->distance_between_neighboring_LHs(26, 27);
                    this->distance_between_neighboring_LHs(25, 27);
                    this->distance_between_neighboring_LHs(24, 27);


                    this->distance_within_single_LH(27, 27);
                    this->distance_within_single_LH(24, 27);
                    this->distance_within_single_LH(21, 27);
                    this->distance_within_single_LH(16, 27);
                    this->distance_within_single_LH(11, 27);
                    this->distance_within_single_LH(5, 27);

                    this->distance_within_single_LH(5, 16);
                    this->distance_within_single_LH(5, 21);


          
                    this->distance_between_neighboring_LHs(6, 27);
                    this->distance_between_neighboring_LHs(6, 21);
                    this->distance_between_neighboring_LHs(6, 16);

                    this->distance_within_single_LH(6, 27);                    
                    this->distance_within_single_LH(6, 21);
                    this->distance_within_single_LH(6, 16);
                    this->distance_within_single_LH(6, 11);
                    
                    


                    this->distance_between_neighboring_LHs(5, 5);
                    this->distance_between_neighboring_LHs(6, 6);
                    this->distance_between_neighboring_LHs(11, 11);
                    this->distance_between_neighboring_LHs(16, 16);
                    //this->distance_between_neighboring_LHs(21, 21);

#endif


#ifdef CALCULATE_EVERYTHING
                    this->calculate_gyration_tensor();

                    this->calculate_distance_6_16();
                    this->calculate_distance_6_27();
                    this->calculate_distance_16_27();

                    this->calculate_distance_9_23();
                    this->calculate_distance_9_27();
                    this->calculate_distance_23_27();
#endif

                    



                    vector< row >().swap(this->LH_c);             // erase LH_c
                    vector< row >().swap(this->nucleosome_c);     // erase nucleosome_c
                    vector< row_orientation >().swap(this->nucleosome_orientations);
                    vector< row >().swap(this->DNA_c);            // erase DNA_c
                    vector< row_orientation >().swap(this->DNA_orientations); // erase DNA_orientations
                    vector< row >().swap(this->tail_c);           // erase tail_c
                }

                if ((frame_number % 1000) == 0)
                if (frame_number <= this->starting_frame)
                    cout << ".";
                else
                    cout << "o";
            }
        }
    }

    fclose(simulation_file);
    cout << endl << endl;

    return 0;
}





int trajectory::set_parameters(int cores, int DNA_beads, int LH_beads, 
    int simulated_LH_beads, int tail_beads, int starting_frame, 
    bool save_data, vector < int > &LH_presence, bool symm_linkers){
    this->cores = cores;
    this->DNA_beads = DNA_beads;
    this->LH_beads = LH_beads;
    this->simulated_LH_beads = simulated_LH_beads;
    
    this->tail_beads = tail_beads;

    this->symmetric_linkers = symm_linkers;

    if(this->symmetric_linkers)
        this->lines_per_frame = cores * 4 + DNA_beads * 4 * (cores + 1) + cores*tail_beads + cores*LH_beads;
    else
        this->lines_per_frame = cores * 4 + DNA_beads * 4 * cores + cores*tail_beads + cores*LH_beads;


    this->starting_line_for_nucleosome = 1;
    this->ending_line_for_nucleosome = cores * 4;

    this->starting_line_for_DNA = cores * 4 + 1;
    if(this->symmetric_linkers){

        this->ending_line_for_DNA = cores * 4 + (cores + 1) * DNA_beads * 4;

        this->starting_line_for_tails = cores * 4 + (cores + 1) * DNA_beads * 4 + 1;
        this->ending_line_for_tails = cores * 4 + (cores + 1) * DNA_beads * 4 + cores*tail_beads;

        this->starting_line_for_LH = cores * 4 + DNA_beads * 4 * (cores + 1) + cores*tail_beads + 1;
        this->ending_line_for_LH = cores * 4 + DNA_beads * 4 * (cores + 1) + cores*tail_beads + cores*LH_beads;
    } else {

        this->ending_line_for_DNA = cores * 4 + cores * DNA_beads * 4;

        this->starting_line_for_tails = cores * 4 + cores * DNA_beads * 4 + 1;
        this->ending_line_for_tails = cores * 4 + cores * DNA_beads * 4 + cores*tail_beads;

        this->starting_line_for_LH = cores * 4 + DNA_beads * 4 * cores + cores*tail_beads + 1;
        this->ending_line_for_LH = cores * 4 + DNA_beads * 4 * cores + cores*tail_beads + cores*LH_beads;

    }


    this->starting_frame = starting_frame;
    this->save_data = save_data;

    if (LH_presence.size() > 0){
        this->LH_presence.swap(LH_presence);

    }
    else{
        for (int i = 0; i < cores; i++)
            this->LH_presence.push_back(1);
    }

    return 0;
}



int trajectory::save_to_file(vector < double > values, string name){
    double temp;
    ofstream outfile;
    string file_name = name;

    
    temp = arithmetic_mean(values);
    //cout << "LH-LH contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_LH_interactions, temp) << endl;
    outfile.open(file_name.append(".txt").c_str(), std::ios_base::app);
    outfile << fixed << showpoint;
    outfile.precision(2);
    outfile << temp << endl;
    outfile.close();
    file_name = name;
    outfile.open(file_name.append("_std").append(".txt").c_str(), std::ios_base::app);
    outfile << fixed << showpoint;
    outfile.precision(2);
    outfile << standard_deviation(values, temp) << endl;
    outfile.close();


    
    return 0;
}


int trajectory::calculate_LH_LH_interactions(){
    int i, j, k, l;
    int core_counter_i = 0, core_counter_k = 0;
    double radius_i = 0, radius_k = 0;
    double dr, dx, dy, dz;
    bool skip_LH_bead = false;
    int how_many_LH_beads_j, how_many_LH_beads_l;
    int count_LH_LH_contacts;

    vector < int > LH_contacts(this->LH_beads - 6, 0);

    for (i = 0; i < this->cores - 1; i++){  // analyzed linker (belonging to core i)

        if (this->LH_presence[i]>0){

            how_many_LH_beads_j = this->simulated_LH_beads;

            count_LH_LH_contacts = 0;
            vector< int >().swap(LH_contacts);
            LH_contacts.resize(this->simulated_LH_beads - 6, 0);
            //LH_contacts.resize(this->LH_beads - 6, 0);

            

            for (j = 6; j < how_many_LH_beads_j; j++){ // its beads


                skip_LH_bead = false;
                core_counter_i = i * this->LH_beads + j;
                radius_i = this->LH_ev_hh[j];

                for (k = i + 1; k < this->cores; k++){

                    if (this->LH_presence[k]>0){

                        if (k != i){  // left if I want to change something (do not want to do reformating)

                            how_many_LH_beads_l = this->simulated_LH_beads;
                            for (l = 6; l < how_many_LH_beads_l; l++){
                                core_counter_k = k*this->LH_beads + l;
                                radius_k = this->LH_ev_hh[l];
                                dx = LH_c[core_counter_i].x - LH_c[core_counter_k].x;
                                dy = LH_c[core_counter_i].y - LH_c[core_counter_k].y;
                                dz = LH_c[core_counter_i].z - LH_c[core_counter_k].z;

                                dr = sqrt(dx*dx + dy*dy + dz*dz);

                                if (dr <= ((radius_i + radius_k) / 2)){
                                    LH_contacts[j - 6] = 1;
                                    skip_LH_bead = true;
                                    break;
                                }

                            } // for (l = 0; l < this->LH_beads; l++){
                        } // if (k != i){          
                        if (skip_LH_bead)
                            break;
                    }
                } // for (k = i+1; k < this->cores; k++){    

            }  //  for (j = 0; j < this->LH_beads; j++){ // its beads


            this->frequency_of_LH_LH_interactions.push_back(arithmetic_mean(LH_contacts));
        }

    }  // for (i = 0; i < this->cores-1; i++){  // analyzed linker (belonging to core i)

    //vector< int >().swap(LH_contacts);

    return 0;
}


int trajectory::calculate_LH_DNA_interactions(){
    int i, j, k, l;
    int core_counter = 0, DNA_counter = 0;
    int how_many_LH_beads;
    double radius_i = 0, radius_k = 0;
    double dr, dx, dy, dz;
    bool skip_LH = false;
    bool interacts_with_parental;
    bool interacts_with_nonparental;
    int start_c;
    int start_l, end_l;
    int start_k, end_k;
    
    vector < int > LH_parental_DNA_contacts(this->LH_beads - 6, 0);
    vector < int > LH_parental_DNA_contacts_first(this->LH_beads - 6, 0);
    vector < int > LH_parental_DNA_contacts_second(this->LH_beads - 6, 0);
    vector < int > LH_nonparental_DNA_contacts(this->LH_beads - 6, 0);


    if (this->symmetric_linkers)
        start_c = 0;
    else
        start_c = 1;

    for (i = start_c; i < this->cores; i++){  // analyzed linker histone belonging to core i

        if (this->LH_presence[i]>0){
        

            LH_parental_DNA_contacts.resize(this->simulated_LH_beads - 6, 0);
            LH_parental_DNA_contacts_first.resize(this->simulated_LH_beads - 6, 0);
            LH_parental_DNA_contacts_second.resize(this->simulated_LH_beads - 6, 0);
            LH_nonparental_DNA_contacts.resize(this->simulated_LH_beads - 6, 0);


            how_many_LH_beads = this->simulated_LH_beads;

            for (j = 6; j < how_many_LH_beads; j++){ // its beads

                interacts_with_parental = false;
                interacts_with_nonparental = false;

                core_counter = i * this->LH_beads + j;
                radius_i = this->LH_ev_hl[j];

                if (this->symmetric_linkers){
                    start_k = 0;
                    end_k = this->cores + 1;
                }else{
                    start_k = 0;
                    end_k = this->cores;
                }

                for (k = start_k; k < end_k; k++){
                    for (l = (k)*this->DNA_beads; l < (k + 1)*this->DNA_beads; l++){
                        DNA_counter = l;
                        radius_k = DNA_ev;

                        dx = LH_c[core_counter].x - DNA_c[DNA_counter].x;
                        dy = LH_c[core_counter].y - DNA_c[DNA_counter].y;
                        dz = LH_c[core_counter].z - DNA_c[DNA_counter].z;

                        dr = sqrt(dx*dx + dy*dy + dz*dz);
                        

                        //if (dr <= ((radius_i + radius_k) / 2)){
                        // division by 4 is introduced to make LH-DNA interactions softer
                        if (dr <= ((radius_i + radius_k) / 2)){
                            if(this->symmetric_linkers){
                                if ((i == k) || ((i + 1) == k)){
                                    if (!interacts_with_parental){
                                        LH_parental_DNA_contacts[j - 6] = 1;
                                        //interacts_with_parental = true;

                                        if (i == k){
                                            LH_parental_DNA_contacts_first[j - 6] = 1;
                                        }
                                        else{
                                            LH_parental_DNA_contacts_second[j - 6] = 1;
                                        }
                                    }

                                }
                                else{
                                    LH_nonparental_DNA_contacts[j - 6] = 1;
                                    interacts_with_nonparental = true;
                                } // if ((i == k) || (i == k + 1)){
                            }else{
                                if (((i-1) == k) || (i == k)){
                                    if (!interacts_with_parental){
                                        LH_parental_DNA_contacts[j - 6] = 1;
                                        //interacts_with_parental = true;

                                        if ((i-1) == k){
                                            LH_parental_DNA_contacts_first[j - 6] = 1;
                                        }
                                        else{
                                            LH_parental_DNA_contacts_second[j - 6] = 1;
                                        }
                                    }

                                }
                                else{
                                    LH_nonparental_DNA_contacts[j - 6] = 1;
                                    interacts_with_nonparental = true;
                                } // if ((i == (k-1)) || (i == k)) {
                            }
                        } // if (dr <= (radius_i + radius_k)){

                        /*if (interacts_with_parental & interacts_with_nonparental)
                            break;*/

                    }  // for (l = (k)*this->DNA_beads; l < (k)*this->DNA_beads + this->DNA_beads; l++){

           
                } // for (k = 0; k < (this->cores+1); k++){            

            }  // for (j = 0; j < this->LH_beads; j++){ // its beads

            this->frequency_of_LH_parental_DNA_interactions.push_back(arithmetic_mean(LH_parental_DNA_contacts));
            this->frequency_of_LH_nonparental_DNA_interactions.push_back(arithmetic_mean(LH_nonparental_DNA_contacts));

            this->frequency_of_LH_parental_DNA_interactions_first.push_back(arithmetic_mean(LH_parental_DNA_contacts_first));
            this->frequency_of_LH_parental_DNA_interactions_second.push_back(arithmetic_mean(LH_parental_DNA_contacts_second));

        }
    }  //  for (i = 0; i < this->cores; i++){  // analyzed linker histone belonging to core i


    return 0;
}


#define DO_NOT_SAVE_CORE_COORDINATES
int trajectory::calculate_core_coordinates(int core_number){
    int i, j;

#ifdef SAVE_CORE_COORDINATES
    ofstream outputfile("core_coordinates_test.txt");
    outputfile << fixed << showpoint;
    outputfile.precision(2);
#endif

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

        outputfile << local_core_particles.coordinates[i][0];
        outputfile << " ";
        outputfile << local_core_particles.coordinates[i][1];
        outputfile << " ";
        outputfile << local_core_particles.coordinates[i][2];
        outputfile << endl;

    }

    outputfile.close();


    return 0;
}



int trajectory::calculate_LH_core_interactions(){
    int i, j, k, l;
    int LH_counter = 0, DNA_counter = 0;
    int how_many_LH_beads;
    double dr, dx, dy, dz;
    double core_distance_x, core_distance_y, core_distance_z;
    double local_LH_bead_x, local_LH_bead_y, local_LH_bead_z;
    double radius_i, radius_k;
    bool interacts_with_parental;
    bool interacts_with_nonparental;

    vector< core > core_vector;

    vector < int > LH_parental_core_contacts(this->simulated_LH_beads - 6, 0);
    vector < int > LH_nonparental_core_contacts(this->simulated_LH_beads - 6, 0);

     for (i = 0; i < this->cores; i++){
        this->calculate_core_coordinates(i);
        core_vector.push_back(local_core_particles);
    }


    for (i = 0; i < this->cores; i++){  // analyzed linker histone  

        if (this->LH_presence[i]>0){

            LH_parental_core_contacts.resize(this->LH_beads - 6, 0);
            LH_nonparental_core_contacts.resize(this->LH_beads - 6, 0);

            how_many_LH_beads = this->simulated_LH_beads;

            for (j = 6; j < how_many_LH_beads; j++){ // the beads of analyzed linker histone

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
                    for (l = 0; l < 300; l++){

                        radius_k = core_ev;
                        dx = local_LH_bead_x - core_vector[k].coordinates[l][0] + core_distance_x;
                        dy = local_LH_bead_y - core_vector[k].coordinates[l][1] + core_distance_y;
                        dz = local_LH_bead_z - core_vector[k].coordinates[l][2] + core_distance_z;
                        dr = sqrt(dx*dx + dy*dy + dz*dz);

                        if (dr <= ((radius_i + radius_k) / 2)){
                            if (i == k){
                                LH_parental_core_contacts[j - 6] = 1;
                                interacts_with_parental = true;
                            }
                            else{
                                LH_nonparental_core_contacts[j - 6] = 1;
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
        

    }

      return 0;
}


int trajectory::calculate_core_core_distance(){
    int i;
    double distance = 0;
    double dist_x, dist_y, dist_z;
    ostringstream file_name;
    ofstream outfile;
    file_name << "core_distances_c" << this->cores << "_b" << this->DNA_beads << ".txt";
    
    outfile.open(file_name.str().c_str(), std::ios_base::app);
    outfile << fixed << showpoint;
    outfile << setw(8);
    outfile.precision(3);

    for (i = 0; i < (this->cores - 1); i++){
        dist_x = this->nucleosome_c[i].x - this->nucleosome_c[i+1].x;
        dist_y = this->nucleosome_c[i].y - this->nucleosome_c[i + 1].y;
        dist_z = this->nucleosome_c[i].z - this->nucleosome_c[i + 1].z;
        distance = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
        outfile << setw(12) << distance << endl;
    }
    outfile.close();
    return 0;
}



// STOPPED HERE - June 6th 2018
int trajectory::calculate_angle_between_cores(){

    int i, index;
    double normal_x = 0;
    double normal_y = 0;
    double normal_z = 1;
    double first_x, first_y, first_z;
    double second_x, second_y, second_z;
    double angle;    
    row first, second;

    ostringstream file_name;
    ofstream outfile;
    file_name << "angle_between_cores_c" << this->cores << "_b" << this->DNA_beads << ".txt";

    outfile.open(file_name.str().c_str(), std::ios_base::app);
    outfile << fixed << showpoint;
    outfile << setw(8);
    outfile.precision(3);

    for (i = 0; i < (this->cores - 1); i++){
        
        index = i;
        first.x = nucleosome_orientations[index].a*normal_x +
            nucleosome_orientations[index + 1].a*normal_y +
            nucleosome_orientations[index + 2].a*normal_z;

        first.y = nucleosome_orientations[index].b*normal_x +
            nucleosome_orientations[index + 1].b*normal_y +
            nucleosome_orientations[index + 2].b*normal_z;

        first.z = nucleosome_orientations[index].c*normal_x +
            nucleosome_orientations[index + 1].c*normal_y +
            nucleosome_orientations[index + 2].c*normal_z;


        index = i + 3;
        second.x = nucleosome_orientations[index].a*normal_x +
            nucleosome_orientations[index + 1].a*normal_y +
            nucleosome_orientations[index + 2].a*normal_z;

        second.y = nucleosome_orientations[index].b*normal_x +
            nucleosome_orientations[index + 1].b*normal_y +
            nucleosome_orientations[index + 2].b*normal_z;

        second.z = nucleosome_orientations[index].c*normal_x +
            nucleosome_orientations[index + 1].c*normal_y +
            nucleosome_orientations[index + 2].c*normal_z;

        angle = 180 * atan2(normr(cross_product(first, second)), dot_product(first, second)) / PI;

        this->angle_between_cores.push_back(angle);

        outfile << setw(12) << angle << endl;

    }
    outfile.close();

    return 0;


}


# define CALC_DIST
int trajectory::calculate_DNA_bending_angles(){

    int i, first_index, second_index;
    double x_first_start, y_first_start, z_first_start;
    double x_first_end, y_first_end, z_first_end;

    double x_second_start, y_second_start, z_second_start;
    double x_second_end, y_second_end, z_second_end;
    double X_1_local, Y_1_local, Z_1_local;
    double X_end_local, Y_end_local, Z_end_local;
    double dist;
    double temp_angle, angle;
    int end_core_i;
    row V21, V34, V43;

    ostringstream file_name;
    ofstream outfile;
    file_name << "DNA_bending_angles_c" << this->cores << "_b" << this->DNA_beads << ".txt";

    outfile.open(file_name.str().c_str(), std::ios_base::app);
    outfile << fixed << showpoint;
    outfile << setw(8);
    outfile.precision(3);

    if (this->symmetric_linkers)
        end_core_i = this->cores;
    else
        end_core_i = this->cores-1;

    for (i = 0; i <= end_core_i; i++){  // analyzed DNA linker an its beads

        if (i == 0){
            first_index = i*this->DNA_beads;
            x_first_start = this->DNA_c[first_index].x;
            y_first_start = this->DNA_c[first_index].y;
            z_first_start = this->DNA_c[first_index].z;

            x_first_end = this->DNA_c[first_index + 1].x;
            y_first_end = this->DNA_c[first_index + 1].y;
            z_first_end = this->DNA_c[first_index + 1].z;
#ifdef CALC_DIST
            dist = sqrt((x_first_start - x_first_end)*(x_first_start - x_first_end) +
                (y_first_start - y_first_end)*(y_first_start - y_first_end) +
                (z_first_start - z_first_end)*(z_first_start - z_first_end));
#endif
        }
        else{
            X_end_local = nucleosome_orientations[(i - 1)*3].a*X_end +
                nucleosome_orientations[(i - 1)*3 + 1].a*Y_end +
                nucleosome_orientations[(i - 1)*3 + 2].a*Z_end;

            Y_end_local = nucleosome_orientations[(i - 1)*3].b*X_end +
                nucleosome_orientations[(i - 1)*3 + 1].b*Y_end +
                nucleosome_orientations[(i - 1)*3 + 2].b*Z_end;

            Z_end_local = nucleosome_orientations[(i - 1)*3].c*X_end +
                nucleosome_orientations[(i - 1)*3 + 1].c*Y_end +
                nucleosome_orientations[(i - 1)*3 + 2].c*Z_end;

            x_first_start = X_end_local + this->nucleosome_c[i-1].x;
            y_first_start = Y_end_local + this->nucleosome_c[i-1].y;
            z_first_start = Z_end_local + this->nucleosome_c[i-1].z;            

            first_index = i*this->DNA_beads; 
            x_first_end = this->DNA_c[first_index].x;
            y_first_end = this->DNA_c[first_index].y;
            z_first_end = this->DNA_c[first_index].z;

#ifdef CALC_DIST
            dist = sqrt((x_first_start - x_first_end)*(x_first_start - x_first_end) +
                (y_first_start - y_first_end)*(y_first_start - y_first_end) +
                (z_first_start - z_first_end)*(z_first_start - z_first_end));
#endif
        }



        if (i == this->cores){

            second_index = (i + 1)*this->DNA_beads - 1;
            x_second_start = this->DNA_c[second_index - 1].x;
            y_second_start = this->DNA_c[second_index - 1].y;
            z_second_start = this->DNA_c[second_index - 1].z;

            x_second_end = this->DNA_c[second_index].x;
            y_second_end = this->DNA_c[second_index].y;
            z_second_end = this->DNA_c[second_index].z;
#ifdef CALC_DIST
            dist = sqrt((x_second_start - x_second_end)*(x_second_start - x_second_end) +
                (y_second_start - y_second_end)*(y_second_start - y_second_end) +
                (z_second_start - z_second_end)*(z_second_start - z_second_end));
#endif
        }
        else{
            second_index = (i + 1)*this->DNA_beads - 1;
            
            x_second_start = this->DNA_c[second_index].x;
            y_second_start = this->DNA_c[second_index].y;
            z_second_start = this->DNA_c[second_index].z;

        
            X_1_local = nucleosome_orientations[i * 3].a*X_1 +
                nucleosome_orientations[i * 3 + 1].a*Y_1 +
                nucleosome_orientations[i * 3 + 2].a*Z_1;

            Y_1_local = nucleosome_orientations[i * 3].b*X_1 +
                nucleosome_orientations[i * 3 + 1].b*Y_1 +
                nucleosome_orientations[i * 3 + 2].b*Z_1;

            Z_1_local = nucleosome_orientations[i * 3].c*X_1 +
                nucleosome_orientations[i * 3 + 1].c*Y_1 +
                nucleosome_orientations[i * 3 + 2].c*Z_1;

            x_second_end = X_1_local + this->nucleosome_c[i].x;
            y_second_end = Y_1_local + this->nucleosome_c[i].y;
            z_second_end = Z_1_local + this->nucleosome_c[i].z;
#ifdef CALC_DIST
            dist = sqrt((x_second_start - x_second_end)*(x_second_start - x_second_end) +
                (y_second_start - y_second_end)*(y_second_start - y_second_end) +
                (z_second_start - z_second_end)*(z_second_start - z_second_end));
#endif

        }

        V21.x = x_first_end - x_first_start;
        V21.y = y_first_end - y_first_start;
        V21.z = z_first_end - z_first_start;

        V34.x = x_second_start - x_second_end;
        V43.x = -V34.x;
        V34.y = y_second_start - y_second_end;
        V43.y = -V34.y;
        V34.z = z_second_start - z_second_end;
        V43.z = -V34.z;

       
        angle = 180 - 180 * atan2(normr(cross_product(V21, V43)), dot_product(V21, V43)) / PI;

        this->bending_angle.push_back(angle);

        outfile << setw(12) << angle << endl;
    }

    outfile.close();

    return 0;
}



int trajectory::calculate_angle_between_DNAs(){

    int i, first_index, second_index;
    row first_DNA_start, first_DNA_end, second_DNA_start, second_DNA_end;
    row first_DNA, second_DNA;
    //double product, norm;
    double angle;
    int start_c;

    char buffer[256];
    string file_name = "entry_exit_angle_cores_";
    sprintf(buffer, "%d_beads_%d.txt", this->cores, this->DNA_beads);
    file_name.append(buffer);
    fstream outfile;

    outfile.open(file_name.c_str(), fstream::app | fstream::out);
    outfile << fixed << showpoint;
    outfile.precision(2);

    if (this->symmetric_linkers)
        start_c = 0;
    else
        start_c = 1;

    for (i = start_c; i < this->cores; i++){

       
        if(this->symmetric_linkers)
            first_index = i*this->DNA_beads;
        else
            first_index = (i-1)*this->DNA_beads;

        first_DNA_start.x = this->DNA_c[first_index].x;
        first_DNA_start.y = this->DNA_c[first_index].y;
        first_DNA_start.z = this->DNA_c[first_index].z;

        if (this->symmetric_linkers)
            first_index = (i + 1)*this->DNA_beads - 1;
        else
         first_index = (i)*this->DNA_beads - 1;
         first_DNA_end.x = this->DNA_c[first_index].x;
         first_DNA_end.y = this->DNA_c[first_index].y;
         first_DNA_end.z = this->DNA_c[first_index].z;

        if (this->symmetric_linkers)
            second_index = (i + 1)*this->DNA_beads;
        else
            second_index = (i)*this->DNA_beads;
        second_DNA_start.x = this->DNA_c[second_index].x;
        second_DNA_start.y = this->DNA_c[second_index].y;
        second_DNA_start.z = this->DNA_c[second_index].z;

       
        if (this->symmetric_linkers)
            second_index = (i + 2)*this->DNA_beads - 1;
        else
            second_index = (i + 1)*this->DNA_beads - 1;

        second_DNA_end.x = this->DNA_c[second_index].x;
        second_DNA_end.y = this->DNA_c[second_index].y;
        second_DNA_end.z = this->DNA_c[second_index].z;

        first_DNA.x = first_DNA_start.x - first_DNA_end.x;
        first_DNA.y = first_DNA_start.y - first_DNA_end.y;
        first_DNA.z = first_DNA_start.z - first_DNA_end.z;

        second_DNA.x = second_DNA_end.x - second_DNA_start.x;
        second_DNA.y = second_DNA_end.y - second_DNA_start.y;
        second_DNA.z = second_DNA_end.z - second_DNA_start.z;


        /*product = dot_product(first_DNA, second_DNA);
        norm = sqrt(first_DNA.x*first_DNA.x + first_DNA.y*first_DNA.y + first_DNA.z*first_DNA.z);
        norm = norm*sqrt(second_DNA.x*second_DNA.x + second_DNA.y*second_DNA.y + second_DNA.z*second_DNA.z);

        product = product / norm;
        product = 180 * acos(product) / PI;*/

        angle = 180 * atan2(normr(cross_product(first_DNA, second_DNA)), dot_product(first_DNA, second_DNA)) / PI;

        outfile << angle << endl;
        this->angle_between_DNAs.push_back(angle);
        

       }

    outfile.close();


    return 0;
}

int trajectory::calculate_1_11_22_angle(){

    row first_vector, second_vector;
    int i;
    int LH_counter = 0;
    int how_many_LH_beads;
    double product, norm;
    double angle;

    char buffer[256];
    string file_name = "angle_1_11_22_cores_";
    sprintf(buffer, "%d_beads_%d.txt", this->cores, this->DNA_beads);
    file_name.append(buffer);
    fstream outfile;

    outfile.open(file_name.c_str(), fstream::app | fstream::out);
    outfile << fixed << showpoint;
    outfile << setw(8);
    outfile.precision(2);
        

    for (i = 0; i < this->cores; i++){  // analyzed linker histone

        how_many_LH_beads = this->simulated_LH_beads;
        if (this->LH_presence[i]>0){
            //first_vector.x = this->
            LH_counter = i * 28;

            first_vector.x = LH_c[LH_counter + 6].x - LH_c[LH_counter + 16].x;
            first_vector.y = LH_c[LH_counter + 6].y - LH_c[LH_counter + 16].y;
            first_vector.z = LH_c[LH_counter + 6].z - LH_c[LH_counter + 16].z;

            second_vector.x = LH_c[LH_counter + how_many_LH_beads - 1].x - LH_c[LH_counter + 16].x;
            second_vector.y = LH_c[LH_counter + how_many_LH_beads - 1].y - LH_c[LH_counter + 16].y;
            second_vector.z = LH_c[LH_counter + how_many_LH_beads - 1].z - LH_c[LH_counter + 16].z;

       

            angle = 180 * atan2(normr(cross_product(first_vector, second_vector)), dot_product(first_vector, second_vector)) / PI;

            if (isnan(angle) == 0){
                this->angle_between_1_11_22.push_back(angle);
                outfile << setw(8) << angle << endl;
                //cout << angle << endl;
            }
            else{
                cout << "Nan" << endl;
                outfile << setw(8) << 1000.0 << endl;
            }
        }
    }

    outfile.close();

    return 0;
}


int trajectory::distance_within_single_LH(int b1, int b2){
    int i;
    int core_counter_1, core_counter_2;
    vector < double > distance;
    double dr, dx, dy, dz;
    char buffer[256];
    string file_name = "within_LH_cores_";
    sprintf(buffer, "%d_beads_%d_distance_%d_%d.txt", this->cores, this->DNA_beads, b1, b2);
    file_name.append(buffer);
    ofstream outfile;
    core_counter_1 = 0;
    core_counter_2 = 0;

    for (i = 0; i < this->cores; i++){
        if (this->LH_presence[i]>0){
            core_counter_1 = i * this->LH_beads + b1;
            core_counter_2 = i * this->LH_beads + b2;
            dx = LH_c[core_counter_1].x - LH_c[core_counter_2].x;
            dy = LH_c[core_counter_1].y - LH_c[core_counter_2].y;
            dz = LH_c[core_counter_1].z - LH_c[core_counter_2].z;
            dr = sqrt(dx*dx + dy*dy + dz*dz);
            distance.push_back(dr);
        }
    }

    outfile.open(file_name.c_str(), std::ios_base::app);
    outfile << fixed << showpoint;
    outfile.precision(2);
    for (vector< double>::iterator it = distance.begin(); it != distance.end(); ++it){
        outfile << *it << endl;
    }
    outfile.close();
    vector< double >().swap(distance);

    return 0;
}


// STOPPED HERE

int trajectory::distance_between_neighboring_LHs(int b1, int b2){

    int i, j;
    int core_counter_i = 0, core_counter_j = 0;
    vector < double > distance;
    double dr, dx, dy, dz;
    char buffer[256];
    string file_name = "cores_";
    sprintf(buffer, "%d_beads_%d_distance_%d_%d.txt", this->cores, this->DNA_beads, b1, b2);
    file_name.append(buffer);
    ofstream outfile;

    for (i = 0; i < (this->cores-1); i++){

        if (this->LH_presence[i]>0){

            core_counter_i = i * this->LH_beads + b1;
            for (j = i+1; j < this->cores; j++){
                if (this->LH_presence[j]>0){

                    if (i != j){
                        core_counter_j = j * this->LH_beads + b2;
                        dx = LH_c[core_counter_i].x - LH_c[core_counter_j].x;
                        dy = LH_c[core_counter_i].y - LH_c[core_counter_j].y;
                        dz = LH_c[core_counter_i].z - LH_c[core_counter_j].z;
                        dr = sqrt(dx*dx + dy*dy + dz*dz);
                        distance.push_back(dr);
                    }
                }
            }
        }
    }

    outfile.open(file_name.c_str(), std::ios_base::app);
    outfile << fixed << showpoint;
    outfile.precision(2);    

    for (vector< double>::iterator it = distance.begin(); it != distance.end(); ++it){
        outfile << *it << endl;
    }
    outfile.close();
    vector< double >().swap(distance);

    return 0;
}

int trajectory::radius_of_gyration_per_core(){

    int i, j;
    int number_LHs = 0;

    for (std::vector<int>::iterator it = this->LH_presence.begin(); it != this->LH_presence.end(); ++it){
        if (*it > 0)
            number_LHs++;
    }
    
          

    std::vector<double> radius_of_gyration_per_core(number_LHs);
    std::vector<double> radius_of_gyration_per_core_std(number_LHs);

    string file_name;
    ofstream outfile;
    char numstr[50];

    
    for (i = 0; i < this->rg.size(); i = i + number_LHs){
        for (j = 0; j < number_LHs; j++){
            radius_of_gyration_per_core[j] = radius_of_gyration_per_core[j] + this->rg[i + j];
        }
    }


   
    for (j = 0; j < number_LHs; j++){
        radius_of_gyration_per_core[j] = radius_of_gyration_per_core[j] / (this->rg.size() / number_LHs);
        file_name = "RG_cores_";        
        //sprintf(numstr, "%d", this->cores);
        sprintf(numstr, "%d", number_LHs);
        file_name = file_name + numstr;
        sprintf(numstr, "%d", j);
        file_name = file_name + string("_core_") + numstr + string(".txt");
        outfile.open(file_name.c_str(), std::ios_base::app);
        outfile << fixed << showpoint;
        outfile.precision(2);
        outfile << radius_of_gyration_per_core[j] << endl;
        outfile.close();
    }

   
    for (i = 0; i < this->rg.size(); i = i + number_LHs){
        for (j = 0; j < number_LHs; j++){
            radius_of_gyration_per_core_std[j] = radius_of_gyration_per_core_std[j] + 
                (this->rg[i + j] - radius_of_gyration_per_core[j])*
                (this->rg[i + j] - radius_of_gyration_per_core[j]);
        }
    }

   
    for (j = 0; j < number_LHs; j++){
        radius_of_gyration_per_core_std[j] = sqrt(radius_of_gyration_per_core_std[j] / ((this->rg.size() / number_LHs) - 1));
        file_name = "RG_cores_";        
        sprintf(numstr, "%d", number_LHs);
        file_name = file_name + numstr;
        sprintf(numstr, "%d", j);
        file_name = file_name + string("_core_") + numstr + string("_std.txt");
        outfile.open(file_name.c_str(), std::ios_base::app);
        outfile << fixed << showpoint;
        outfile.precision(2);
        outfile << radius_of_gyration_per_core_std[j] << endl;
        outfile.close();
    }


    return 0;
}






    int trajectory::calculate_LH_tail_interactions(){
        int i, j, k, l;
        int core_counter = 0, tail_counter = 0;
        double radius_i = 0, radius_k = 0;
        double dr, dx, dy, dz;
        bool skip_LH = false;


        vector < int > LH_parental_H2A1_contacts(this->simulated_LH_beads - 6, 0);
        vector < int > LH_nonparental_H2A1_contacts(this->simulated_LH_beads - 6, 0);

        vector < int > LH_parental_H2A2_contacts(this->simulated_LH_beads - 6, 0);
        vector < int > LH_nonparental_H2A2_contacts(this->simulated_LH_beads - 6, 0);

        vector < int > LH_parental_H2B_contacts(this->simulated_LH_beads - 6, 0);
        vector < int > LH_nonparental_H2B_contacts(this->simulated_LH_beads - 6, 0);

        vector < int > LH_parental_H3_contacts(this->simulated_LH_beads - 6, 0);
        vector < int > LH_nonparental_H3_contacts(this->simulated_LH_beads - 6, 0);

        vector < int > LH_parental_H4_contacts(this->simulated_LH_beads - 6, 0);
        vector < int > LH_nonparental_H4_contacts(this->simulated_LH_beads - 6, 0);



        radius_k = CDT_tail;

        for (i = 0; i < this->cores; i++){  // analyzed linker histone belonging to core i

            if (this->LH_presence[i]>0){


                LH_parental_H2A1_contacts.resize(this->simulated_LH_beads - 6, 0);
                LH_nonparental_H2A1_contacts.resize(this->simulated_LH_beads - 6, 0);

                LH_parental_H2A2_contacts.resize(this->simulated_LH_beads - 6, 0);
                LH_nonparental_H2A2_contacts.resize(this->simulated_LH_beads - 6, 0);

                LH_parental_H2B_contacts.resize(this->simulated_LH_beads - 6, 0);
                LH_nonparental_H2B_contacts.resize(this->simulated_LH_beads - 6, 0);

                LH_parental_H3_contacts.resize(this->simulated_LH_beads - 6, 0);
                LH_nonparental_H3_contacts.resize(this->simulated_LH_beads - 6, 0);

                LH_parental_H4_contacts.resize(this->simulated_LH_beads - 6, 0);
                LH_nonparental_H4_contacts.resize(this->simulated_LH_beads - 6, 0);

              

                   for (j = 6; j < this->simulated_LH_beads; j++){ // its beads

                    radius_i = this->LH_ev_ht[j];
                    core_counter = i * this->LH_beads + j;

                    for (k = 0; k < this->cores; k++){
                        for (l = (k)*this->tail_beads; l < (k + 1)*this->tail_beads; l++){

                            tail_counter = l;

                            dx = LH_c[core_counter].x - tail_c[tail_counter].x;
                            dy = LH_c[core_counter].y - tail_c[tail_counter].y;
                            dz = LH_c[core_counter].z - tail_c[tail_counter].z;
                            dr = sqrt(dx*dx + dy*dy + dz*dz);

                            if (dr < ((radius_i + radius_k) / 2)){


                                if (((l % this->tail_beads) > 25) & ((l % this->tail_beads) < 34)){ // H2A1
                                    //cout << l << " H2A1" << endl;

                                    if (i == k){
                                        LH_parental_H2A1_contacts[j - 6] = 1;
                                    }
                                    else{
                                        LH_nonparental_H2A1_contacts[j - 6] = 1;
                                    }
                                } // if (((l % this->tail_beads) < 34) & ((l % this->tail_beads) > 25)){


                                if (((l % this->tail_beads) > 43) & ((l % this->tail_beads) < 50)){ // H2A2
                                    

                                    if (i == k){
                                        LH_parental_H2A2_contacts[j - 6] = 1;
                                    }
                                    else{
                                        LH_nonparental_H2A2_contacts[j - 6] = 1;
                                    }
                           


                                if (((l % this->tail_beads) > 33) & ((l % this->tail_beads) < 44)){ // H2B
                                    
                                    if (i == k){
                                        LH_parental_H2B_contacts[j - 6] = 1;
                                    }
                                    else{
                                        LH_nonparental_H2B_contacts[j - 6] = 1;
                                    }
                              

                                if ((l % this->tail_beads) < 16){ // H3
                                    //cout << l << " H3" << endl;

                                    if (i == k){
                                        LH_parental_H3_contacts[j - 6] = 1;
                                    }
                                    else{
                                        LH_nonparental_H3_contacts[j - 6] = 1;
                                    }
                                } // if ((l % this->tail_beads) < 16){ // H3


                                if (((l % this->tail_beads) > 15) & ((l % this->tail_beads) < 26)){ // H4
                                    //cout << l << " H4" << endl;

                                    if (i == k){
                                        LH_parental_H4_contacts[j - 6] = 1;
                                    }
                                    else{
                                        LH_nonparental_H4_contacts[j - 6] = 1;
                                    }
                        
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
            }


        } // for (i = 0; i < this->cores; i++){  // analyzed linker histone belonging to core i

        return 0;
    }





int trajectory::calculate_statistics(){    
    ofstream outfile;
    double temp;

    cout << fixed << showpoint;
    cout.precision(2);
  
      
    
 
    cout << endl;
    temp = arithmetic_mean(this->rg);
    cout << "radius of gyration mean = " << temp << ", std = " << standard_deviation(this->rg, temp) << endl;
    this->save_to_file(this->rg, "radius_of_gyration");


    cout << endl;
    temp = arithmetic_mean(this->frequency_of_LH_parental_DNA_interactions);
    cout << "Parental LH-DNA contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_parental_DNA_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_parental_DNA_interactions, "LH_parental_DNA");


    temp = arithmetic_mean(this->frequency_of_LH_nonparental_DNA_interactions);
    cout << "Non-parental LH-DNA contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_nonparental_DNA_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_nonparental_DNA_interactions, "LH_nonparental_DNA");

    cout << endl;
    temp = arithmetic_mean(this->frequency_of_LH_parental_DNA_interactions_first);
    cout << "Parental LH-DNA contacts (first) mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_parental_DNA_interactions_first, temp) << endl;
    this->save_to_file(this->frequency_of_LH_parental_DNA_interactions_first, "LH_parental_DNA_first");

    cout << endl;
    temp = arithmetic_mean(this->frequency_of_LH_parental_DNA_interactions_second);
    cout << "Parental LH-DNA contacts (second) mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_parental_DNA_interactions_second, temp) << endl;
    this->save_to_file(this->frequency_of_LH_parental_DNA_interactions_second, "LH_parental_DNA_second");

    
    temp = arithmetic_mean(this->frequency_of_LH_LH_interactions);
    cout << "LH-LH contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_LH_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_LH_interactions, "LH_LH");
    
    
    cout << endl;
    temp = arithmetic_mean(this->frequency_of_LH_parental_core_interactions);
    cout << "Parental LH-core contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_parental_core_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_parental_core_interactions, "LH_parental_core");
    

    temp = arithmetic_mean(this->frequency_of_LH_nonparental_core_interactions);
    cout << "Non-parental LH-core contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_nonparental_core_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_nonparental_core_interactions, "LH_nonparental_core");
    
    //cout << "----------------------------------------------------------------------------------------" << endl;


    cout << endl;
    temp = arithmetic_mean(this->frequency_of_LH_parental_H2A1_interactions);
    cout << "Parental LH-H2A1 contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_parental_H2A1_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_parental_H2A1_interactions, "LH_parental_H2A1");
    

    temp = arithmetic_mean(this->frequency_of_LH_nonparental_H2A1_interactions);
    cout << "Non-parental LH-H2A1 contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_nonparental_H2A1_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_nonparental_H2A1_interactions, "LH_nonparental_H2A1");
    

    cout << endl;
    temp = arithmetic_mean(this->frequency_of_LH_parental_H2A2_interactions);
    cout << "Parental LH-H2A2 contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_parental_H2A2_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_parental_H2A2_interactions, "LH_parental_H2A2");
    


    temp = arithmetic_mean(this->frequency_of_LH_nonparental_H2A2_interactions);
    cout << "Non-parental LH-H2A2 contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_nonparental_H2A2_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_nonparental_H2A2_interactions, "LH_nonparental_H2A2");
    

    cout << endl;
    temp = arithmetic_mean(this->frequency_of_LH_parental_H2B_interactions);
    cout << "Parental LH-H2B contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_parental_H2B_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_parental_H2B_interactions, "LH_parental_H2B");
    


    temp = arithmetic_mean(this->frequency_of_LH_nonparental_H2B_interactions);
    cout << "Non-parental LH-H2B contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_nonparental_H2B_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_nonparental_H2B_interactions, "LH_nonparental_H2B");
    
    cout << endl;
    temp = arithmetic_mean(this->frequency_of_LH_parental_H3_interactions);
    cout << "Parental LH-H3 contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_parental_H3_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_parental_H3_interactions, "LH_parental_H3");
    

    temp = arithmetic_mean(this->frequency_of_LH_nonparental_H3_interactions);
    cout << "Non-parental LH-H3 contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_nonparental_H3_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_nonparental_H3_interactions, "LH_nonparental_H3");
    

    cout << endl;
    temp = arithmetic_mean(this->frequency_of_LH_parental_H4_interactions);
    cout << "Parental LH-H4 contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_parental_H4_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_parental_H4_interactions, "LH_parental_H4");
    

    temp = arithmetic_mean(this->frequency_of_LH_nonparental_H4_interactions);
    cout << "Non-parental LH-H4 contacts mean = " << temp << ", std = " << standard_deviation(this->frequency_of_LH_nonparental_H4_interactions, temp) << endl;
    this->save_to_file(this->frequency_of_LH_nonparental_H4_interactions, "LH_nonparental_H4");
    

    cout << endl;
    temp = arithmetic_mean(this->bending_angle);
    cout << "DNA bending angle mean = " << temp << ", std = " << standard_deviation(this->bending_angle, temp) << endl;      
    this->save_to_file(this->bending_angle, "DNA_bending_angle");

    cout << endl;
    temp = arithmetic_mean(this->angle_between_cores);
    cout << "Angle between cores mean = " << temp << ", std = " << standard_deviation(this->angle_between_cores, temp) << endl;
    this->save_to_file(this->bending_angle, "Angle_between_cores");
    


    temp = arithmetic_mean(this->angle_between_DNAs);
    cout << "Angle between DNAs is " << arithmetic_mean(this->angle_between_DNAs) << ", std = " << standard_deviation(this->angle_between_DNAs, temp) << endl;
    this->save_to_file(this->angle_between_DNAs, "Angle_between_DNAs");
    


    cout << endl;
    temp = arithmetic_mean(this->rg);
    cout << "radius of gyration mean = " << temp << ", std = " << standard_deviation(this->rg, temp) << endl;
    this->save_to_file(this->rg, "radius_of_gyration");

        
    cout << endl;
    temp = arithmetic_mean(this->distance_9_27);
    cout << "distance A = " << temp << ", std = " << standard_deviation(this->distance_9_27, temp) << endl;
    this->save_to_file(this->distance_9_27, "distance_A_9_27");

    cout << endl;
    temp = arithmetic_mean(this->distance_9_23);
    cout << "distance B = " << temp << ", std = " << standard_deviation(this->distance_9_23, temp) << endl;
    this->save_to_file(this->distance_9_23, "distance_B_9_23");

    cout << endl;
    temp = arithmetic_mean(this->distance_23_27);
    cout << "distance C = " << temp << ", std = " << standard_deviation(this->distance_23_27, temp) << endl;
    this->save_to_file(this->distance_23_27, "distance_C_23_27");


    cout << endl;
    temp = arithmetic_mean(this->distance_6_16);
    cout << "1-11 distance mean is  " << temp << ", std = " << standard_deviation(this->distance_6_16, temp) << endl;
    this->save_to_file(this->distance_6_16, "distance_6_16");

    temp = arithmetic_mean(this->distance_6_27);
    cout << "1-22 distance mean is  " << temp << ", std = " << standard_deviation(this->distance_6_27, temp) << endl;
    this->save_to_file(this->distance_6_27, "distance_6_27");

    temp = arithmetic_mean(this->distance_16_27);
    cout << "11-22 distance mean is  " << temp << ", std = " << standard_deviation(this->distance_16_27, temp) << endl;
    this->save_to_file(this->distance_16_27, "distance_16_27");

    radius_of_gyration_per_core();

    return 0;
}
