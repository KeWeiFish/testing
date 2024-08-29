#include "Test Function.h"

void readParameters(const string &filename, vector<int> &Dimension, int &Evaluation, int &Iteration, vector<double> &PopulationRatial, double &F, double &CR, int &functionumber) {
    ifstream inputFile(filename);
    // Check if the file was opened successfully
    if (!inputFile) {
        throw runtime_error("Error opening file: " + filename);
    }

    string line;
    while (getline(inputFile, line)) {
        stringstream ss(line);
        string type;
        ss >> type;

        if (type == "Dimension") {
            int value;
            while (ss >> value) {
                Dimension.push_back(value);
            }
        } 
        else if (type == "Evaluation") 
            ss >> Evaluation;
        else if (type == "Iteration") 
            ss >> Iteration;
        else if (type == "PopulationRatial") {
            int value;
            while (ss >> value) {
                PopulationRatial.push_back(value);
            }
        } 
        else if (type == "F") 
            ss >> F;
        else if (type == "CR") 
            ss >> CR;
        else if (type == "FunctionNumber") {
            ss >> functionumber;
        }
    }
    inputFile.close();
}




void outputResult(const int Dimension, const string FunctionName, const vector<double> IterationGB){
    ostringstream oss;
    oss<<FunctionName<<" in "<<Dimension<<" Dimension"<<endl;
    for(int i=0;i<IterationGB.size();i++)
        oss<<"Iteration "<<i<<": "<<setw(2)<<right<<IterationGB[i]<<endl;
    oss<<string(25, '-')<<endl;
    cout<<oss.str()<<endl;
    return ;
}

map<int, std::string> funcPair = {
    {1,"Ackley"},{2,"Rastrigin"},{3,"HappyCat"},
    {4,"Rosenbrock"},{5,"Zakharov"},{6,"Michalewicz"},
    {7,"Schewfel"},{8,"BentCigar"},{9,"DropWave"},
    {10,"Step"}
    };

void writeResult(const int Dimension, const int funcIndex, const vector<double> data){
    // Construct the filename using the function name and dimension

    __fs::filesystem::path dir{"GBresults"};
    if (!__fs::filesystem::exists(dir))
        __fs::filesystem::create_directory(dir);
    
    __fs::filesystem::path filePath = dir / (funcPair[funcIndex] + "_" + std::to_string(Dimension) + "D.txt");
    // Create an ofstream object for file writing
    ofstream outFile(filePath);

    // Check if the file is open
    if (!outFile.is_open()) {
        cerr << "Failed to open file: " << filePath << endl;
        return;
    }

    for (int i=0;i<data.size();i++) {
        outFile << data[i] << endl;
    }

    outFile.close();
    return ;
}

void writeConvergence(const int Dimension, const int funcIndex, const vector<double> data, const int iter){
    // Construct the filename using the function name and dimension
    __fs::filesystem::path folderpath = funcPair[funcIndex] + "_" + to_string(Dimension) + "D_Convergence";
    __fs::filesystem::path dir{folderpath};

    if (!__fs::filesystem::exists(dir))
        __fs::filesystem::create_directory(dir);
    
    __fs::filesystem::path filePath = dir / (to_string(iter) + "iteraion.txt");
    // Create an ofstream object for file writing
    ofstream outFile(filePath);

    // Check if the file is open
    if (!outFile.is_open()) {
        cerr << "Failed to open file: " << filePath << endl;
        return;
    }

    for (int i=0;i<data.size();i++) {
        outFile << data[i] << endl;
    }

    outFile.close();
    return ;
}