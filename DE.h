
#include "Function.h"

double** create2DArray(const int row, const int column){
    double** temp = new double*[row];
    for (int i = 0; i < row; i++) 
        temp[i] = new double[column];
    return temp;
}

void delete2DArray(double** array, int rows) {
    // Free each row
    for (int i = 0; i < rows; i++) {
        delete[] array[i];
    }
    
    // Free the array of pointers
    delete[] array;
    array = nullptr;
}

void deep_copy(double** a, double** b, const int row, const int column){
    for(int i=0;i<row;i++)
        for(int j=0;j<column;j++)
            a[i][j] = b[i][j];
    return ;
}


double calculate_GB(double** matrix,const int pop, const int dim, const int functionIndex){
    double GB=DBL_MAX,temp=0;
    for(int i=0;i<pop;i++){
        temp = calculate_test_function(matrix[i], dim, functionIndex);
        if(GB>temp)
            GB=temp;
    }
    return GB;
}

struct IndexValue{
    double value;
    int index;
};

bool compare(const IndexValue &a, const IndexValue &b) {
    return a.value > b.value;
}


double calculate_GB_replace_last(double** target, double** first,const int pop, const int dim, const int functionIndex){
    double GB=DBL_MAX,temp=0;
    vector<IndexValue> gblist;
    for(int i=0;i<pop;i++){
        temp = calculate_test_function(target[i], dim, functionIndex);
        gblist.push_back({temp,i});
        if(GB>temp)
            GB=temp;
    }
    sort(gblist.begin(), gblist.end(), compare);

    for(int i=0;i<pop/20;i++)
        for(int j=0;j<dim;j++){
            int test = gblist[i].index;
            target[test][j] = first[test][j]; 
        }

    return GB;
}

double** selection(double** target_matrix, double** trial_matrix, const int pop, const int dim, const int functionIndex){
    double** temp = create2DArray(pop,dim);

    for(int i=0;i<pop;i++){
        if(calculate_test_function(target_matrix[i],dim,functionIndex) <= calculate_test_function(trial_matrix[i],dim,functionIndex))
            for(int j=0;j<dim;j++)
                temp[i][j] = target_matrix[i][j];
        else
            for(int j=0;j<dim;j++)
                temp[i][j] = trial_matrix[i][j];
    }
    return temp;
}


double** crossover(double** target_matrix, double** mutant_matrix, const int pop, const int dim, const double F, const double CR){
    double** trial_matrix = create2DArray(pop,dim);

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> double_dis(0.0, 1.0);
    // Create a uniform distribution between 0 and 1

    int jrand = rand()%dim;
    for(int i=0;i<pop;i++){
        for (int j = 0; j < dim; ++j) {
            if (double_dis(gen)<= CR || j == jrand) {
                trial_matrix[i][j] = mutant_matrix[i][j];
            } else {
                trial_matrix[i][j] = target_matrix[i][j];
            }
        }
    }
    return trial_matrix;
}


double** mutation(double** target_matrix, const int pop, const int dim, const double F, const double GB){
    double** mutant_matrix = create2DArray(pop,dim);

    for(int i=0;i<pop;i++){
        int x1, x2, x3, x4, x5;
        x1 = rand() % pop;
        do {
            x2 = rand() % pop;
        } while (x2 == x1);
                do {
            x3 = rand() % pop;
        } while (x3 == x2 || x3 == x1);
        do {
            x4 = rand() % pop;
        } while (x4 == x3 || x4 == x2 || x4 == x1);
                do {
            x5 = rand() % pop;
        } while (x5 == x4 || x4 == x3 || x4 == x2 || x4 == x1);
        for(int j=0;j<dim;j++){
            mutant_matrix[i][j] = target_matrix[x1][j] + F*(target_matrix[x2][j]-target_matrix[x3][j]) + F*(target_matrix[x4][j]-target_matrix[x5][j]);
        }
    }
    return mutant_matrix;
}


// Function to generate random points in a 2D dynamic array
double** generate_random_points(int collectingpop, int dim, const double minValue, const double maxValue) {
    double** points = new double*[collectingpop];
    random_device rd;
    mt19937 gen(rd());

    for (int i = 0; i < collectingpop; ++i) {
        points[i] = new double[dim];
        for (int j = 0; j < dim; ++j) {
            uniform_real_distribution<> dis(minValue, maxValue);
            points[i][j] = dis(gen);
        }
    }

    return points;
}

// Function to calculate Euclidean distance between two points
double euclidean_distance(const double* a, const double* b, int dim) {
    double distance = 0.0;
    for (int i = 0; i < dim; ++i) {
        distance += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sqrt(distance);
}

// k-means clustering algorithm
void kmeans_Pop(double** matrix, int pop, int dim, int max_iterations, const double minValue, const double maxValue) {
    int collectingpop = pop*100;
    int counts[collectingpop],labels[collectingpop];
    double** points = generate_random_points(collectingpop,dim,minValue,maxValue);
    double** temp_Pop = create2DArray(pop,dim);
    // Initialize matrix by randomly selecting points from points
    random_device rd;
    mt19937 gen(rd());


    // Create a vector containing the range of numbers
    vector<int> numbers(collectingpop - 1);
    std::iota(numbers.begin(), numbers.end(), 0);
    // Shuffle the vector to get a random order
    shuffle(numbers.begin(), numbers.end(), gen);
    // Select the first num_samples elements from the shuffled vector
    vector<int> sampled_numbers(numbers.begin(), numbers.begin() + pop);

    for(int i=0;i<pop;i++)
        for(int j=0;j<dim;j++)
            matrix[i][j] = points[sampled_numbers[i]][j];

    for (int iter = 0; iter < max_iterations; iter++) {
        // Assign points to the nearest center
        for (int i = 0; i < collectingpop; i++) {
            double min_distance = euclidean_distance(points[i], matrix[0], dim);
            int label = 0;
            for (int j = 1; j < pop; j++) {
                double current_distance = euclidean_distance(points[i], matrix[j], dim);
                if (current_distance < min_distance) {
                    min_distance = current_distance;
                    label = j;
                }
            }
            labels[i] = label;
        }

        //Initialize temp_Pop(2D) and counts with value 0
        for (int i = 0; i < pop; ++i) {
            fill(temp_Pop[i], temp_Pop[i] + dim, 0.0);
        }
        fill(counts, counts + collectingpop, 0);

        for (int i = 0; i < collectingpop; ++i) {
            int label = labels[i];
            for (int j = 0; j < dim; ++j) {
                temp_Pop[label][j] += points[i][j];
            }
            counts[label]++;
        }

        for (int i = 0; i < pop; ++i) {
            for (int j = 0; j < dim; ++j) {
                if (counts[i] > 0) {
                    temp_Pop[i][j] /= counts[i];
                }
            }
        }

        // Check for convergence
        double shift = 0.0;
        for (int i = 0; i < pop; ++i) {
            shift += euclidean_distance(matrix[i], temp_Pop[i], dim);
        }

        if (shift < 1e-6) {
            // cout<<iter<<endl;
            break;
        }

        // Update matrix
        for (int i = 0; i < pop; ++i) {
            copy(temp_Pop[i], temp_Pop[i] + dim, matrix[i]);
        }
    }

    // Clean up
    delete2DArray(temp_Pop,pop);
    return ;
}

void PopulationInitialization (double** matrix,const int row, const int column,const double minValue, double const maxValue){
    random_device rd;  // Obtain a random number from hardware
    mt19937 gen(rd()); // Seed the generator
    uniform_real_distribution<> dis(minValue, maxValue); // Define the range

    for (int i = 0; i < row; i++)
        for (int j = 0; j < column; j++)
            matrix[i][j] = dis(gen);

    return ;
}


double DE(int Dim,int Eval, int PopulationRatial, double F,double CR, const int functionIndex, const int iter){
    int Population = Dim*PopulationRatial;
    double upper_bound,lower_bound,GB;
    double **First_Pop = nullptr, **target_matrix = nullptr, **temp = nullptr, **mutation_matrix = nullptr, **trial_matrix = nullptr;

    vector<double> convergence;
    First_Pop = create2DArray(Population,Dim);
    target_matrix = create2DArray(Population,Dim);
    set_search_bound(&upper_bound, &lower_bound, functionIndex); // Assuming func_index + 1 as 1 for example
    //choose standard Population Initialization or kmeans Population Initialization

    // PopulationInitialization(First_Pop, Population, Dim, lower_bound, upper_bound);
    kmeans_Pop(First_Pop,Population,Dim,100,lower_bound,upper_bound);
    deep_copy(target_matrix,First_Pop,Population,Dim);
    GB = calculate_GB(target_matrix,Population,Dim, functionIndex);
    convergence.push_back(GB);

    for(int run=0;run<Eval*Dim;run++){
        mutation_matrix = mutation(target_matrix,Population,Dim,F,GB);
        trial_matrix = crossover(target_matrix,mutation_matrix,Population,Dim,F,CR);
        delete2DArray(mutation_matrix,Population);
        temp = selection(target_matrix,trial_matrix,Population,Dim,functionIndex);
        delete2DArray(trial_matrix,Population);
        deep_copy(target_matrix,temp,Population,Dim);
        delete2DArray(temp,Population);

        //choose salmon-modification or standard GB calculating
        if(run%300==0)
            GB = calculate_GB_replace_last(target_matrix,First_Pop,Population,Dim,functionIndex);
        else
            GB = calculate_GB(target_matrix,Population,Dim, functionIndex);
        // GB = calculate_GB(target_matrix,Population,Dim, functionIndex);
        convergence.push_back(GB);
    }
    writeConvergence(Dim,functionIndex,convergence,iter);
    delete2DArray(target_matrix,Population);
    delete2DArray(First_Pop,Population);
    return GB;
}