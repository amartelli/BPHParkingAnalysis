#include "tensorflow/core/framework/graph.pb.h"
#include "tensorflow/core/framework/tensor.h"

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/framework/tensor.h"
#include "tensorflow/core/lib/io/path.h"

#include "tensorflow/core/graph/default_device.h"

#include <exception>


int main(int argc, char** argv)
{
    if (argc<=1)
    {
        std::cerr<<"\tSyntax: evaluateModel [graph.pb]"<<std::endl;
        return 1;
    }
    
    std::cout<<"Testing model: "<<argv[1]<<std::endl;

    tensorflow::Status status;

    // load it
    tensorflow::GraphDef graphDef;
    status = ReadBinaryProto(tensorflow::Env::Default(), argv[1], &graphDef);
    tensorflow::graph::SetDefaultDevice("/cpu:0", &graphDef);
    
    // check for success
    if (!status.ok())
    {
        throw std::runtime_error("InvalidGraphDef: error while loading graph def: "+status.ToString());
    }
    
    tensorflow::Session* session;
    tensorflow::SessionOptions opts;
    opts.config.set_intra_op_parallelism_threads(1);
    opts.config.set_inter_op_parallelism_threads(1);
    TF_CHECK_OK(tensorflow::NewSession(opts, &session));
    TF_CHECK_OK(session->Create(graphDef));
    
    tensorflow::Tensor feature_tensor(tensorflow::DT_FLOAT, {1,200,24});
    auto features = feature_tensor.tensor<float,3>();
    for (int i = 0; i < 24; ++i)
    {
        features(0,0,i) = 1;
        for (int j = 1; j < 200; ++j)
        {
            features(0,j,i) = 0;
        }
    }
    
    std::vector<tensorflow::Tensor> outputs; 
    TF_CHECK_OK(session->Run(
        {
            {"features",feature_tensor},
            
        }, //input map
        {"prediction"}, //output node names 
        {}, //additional nodes run but not put in outputs
        &outputs
    ));
    for (const auto& tensor: outputs)
    {
        auto tensor_flat = tensor.flat<float>();
        for (int i = 0; i < tensor_flat.size(); ++i)
        {
            std::cout<<i<<": "<<tensor_flat(0)<<std::endl;
        }
    }
    
    return 0;
}

