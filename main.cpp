#include "DE.h"



int main() {
    vector<int> Dimension;
    vector<double> PopulationRatial;
    int Evaluation, Iteration, FunctionNumber;
    double F, CR;

    readParameters("parameter.txt", Dimension, Evaluation, Iteration, PopulationRatial, F, CR,FunctionNumber);

    for(int dim=0;dim<Dimension.size();dim++){
        for(int func=0;func<FunctionNumber;func++){
            vector<double> ObjectiveValue;
            vector<future<double>> futures;
            for(int iter = 0; iter < Iteration; iter++) {
                    // Launch a task asynchronously
                    auto future = async(launch::async, DE, Dimension[dim], Evaluation, PopulationRatial[dim], F, CR, func+1, iter);
                    futures.push_back(std::move(future));
                }
                // Retrieve the results and populate ObjectiveValue
            for(auto& future : futures) {
                ObjectiveValue.push_back(future.get());
            }
            writeResult(Dimension[dim], func+1, ObjectiveValue);
        }
    }

    return 0;
}