// -*- C++ -*-
#ifndef RIVET_RivetONNXrt_HH
#define RIVET_RivetONNXrt_HH

#include <iostream>
#include <functional>
#include <numeric>

#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Tools/Utils.hh"
#include "onnxruntime/onnxruntime_cxx_api.h"


namespace Rivet {

  /// @brief Simple interface class to take care of basic ONNX networks
  ///
  /// See analyses/examples/EXAMPLE_ONNX.cc for how to use this.
  /// @note A node is not a neuron but a single tensor of arbitrary dimension size
  class RivetONNXrt {

  public:

    // Suppress default constructor
    RivetONNXrt() = delete;

    /// Constructor
    RivetONNXrt(const string& filename, const string& runname = "RivetONNXrt") {

      // Set some ORT variables that need to be kept in memory
      _env = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, runname.c_str());

      // Load the model
      Ort::SessionOptions sessionopts;
      _session = std::make_unique<Ort::Session> (*_env, filename.c_str(), sessionopts);

      // Store network hyperparameters (input/output shape, etc.)
      getNetworkInfo();

      MSG_DEBUG(*this);
    }

    /// Given a multi-node input vector, populate and return the multi-node output vector
    vector<vector<float>> compute(vector<vector<float>>& inputs) const {

      /// Check that number of input nodes matches what the model expects
      if (inputs.size() != _inDims.size()) {
        throw("Expected " + to_string(_inDims.size())
              + " input nodes, received " + to_string(inputs.size()));
      }

      // Create input tensor objects from input data
      vector<Ort::Value> ort_input;
      ort_input.reserve(_inDims.size());
      auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
      for (size_t i=0; i < _inDims.size(); ++i) {

        // Check that input data matches expected input node dimension
        if (inputs[i].size() != _inDimsFlat[i]) {
          throw("Expected flattened input node dimension " + to_string(_inDimsFlat[i])
                 + ", received " + to_string(inputs[i].size()));
        }

        ort_input.emplace_back(Ort::Value::CreateTensor<float>(memory_info,
                                                               inputs[i].data(), inputs[i].size(),
                                                               _inDims[i].data(), _inDims[i].size()));
      }

      // retrieve output tensors
      auto ort_output = _session->Run(Ort::RunOptions{nullptr}, _inNames.data(),
                                      ort_input.data(), ort_input.size(),
                                      _outNames.data(), _outNames.size());

      // construct flattened values and return
      vector<vector<float>> outputs; outputs.resize(_outDims.size());
      for (size_t i = 0; i < _outDims.size(); ++i) {
        float* floatarr = ort_output[i].GetTensorMutableData<float>();
        outputs[i].assign(floatarr, floatarr + _outDimsFlat[i]);
      }
      return outputs;
    }

    /// Given a single-node input vector, populate and return the single-node output vector
    vector<float> compute(const vector<float>& inputs) const {
      if (_inDims.size() != 1 || _outDims.size() != 1) {
        throw("This method assumes a single input/output node!");
      }
      vector<vector<float>> wrapped_inputs = { inputs };
      vector<vector<float>> outputs = compute(wrapped_inputs);
      return outputs[0];
    }

    /// Method to check if @a key exists in network metatdata
    bool hasKey(const std::string& key) const {
      Ort::AllocatorWithDefaultOptions allocator;
      return (bool)_metadata->LookupCustomMetadataMapAllocated(key.c_str(), allocator);
    }

    /// Method to retrieve value associated with @a key
    /// from network metadata and return value as type T
    template <typename T,
      typename std::enable_if_t<!is_iterable_v<T> | is_cstring_v<T> >>
    T retrieve(const std::string& key) const {
      Ort::AllocatorWithDefaultOptions allocator;
      Ort::AllocatedStringPtr res = _metadata->LookupCustomMetadataMapAllocated(key.c_str(), allocator);
      if (!res) {
        throw("Key '"+key+"' not found in network metadata!");
      }
      /*if constexpr (std::is_same<T, std::string>::value) {
        return res.get();
      }*/
      return lexical_cast<T>(res.get());
    }

    /// Template specialisation of retrieve for std::string
    std::string retrieve(const std::string& key) const {
      Ort::AllocatorWithDefaultOptions allocator;
      Ort::AllocatedStringPtr res = _metadata->LookupCustomMetadataMapAllocated(key.c_str(), allocator);
      if (!res) {
        throw("Key '"+key+"' not found in network metadata!");
      }
      return res.get();
    }

    /// Overload of retrieve for vector<T>
    template <typename T>
    vector<T> retrieve(const std::string & key) const {
      const vector<string> stringvec = split(retrieve(key), ",");
      vector<T> returnvec = {};
      for (const string & s : stringvec){
        returnvec.push_back(lexical_cast<T>(s));
      }
      return returnvec;
    }

    /// Overload of retrieve for vector<T>, with a default return
    template <typename T>
    vector<T> retrieve(const std::string & key, const vector<T> & defaultreturn) const {
      try {
        return retrieve<T>(key);
      } catch (...) {
        return defaultreturn;
      }
    }

    std::string retrieve(const std::string& key, const std::string& defaultreturn) const {
      try {
        return retrieve(key);
      } catch (...) {
        return defaultreturn;
      }
    }

    /// Variation of retrieve method that falls back
    /// to @a defaultreturn if @a key cannot be found
    template <typename T,
      typename std::enable_if_t<!is_iterable_v<T> | is_cstring_v<T> >>
    T retrieve(const std::string& key, const T& defaultreturn) const {
      try {
        return retrieve<T>(key);
      } catch (...) {
        return defaultreturn;
      }
    }

    /// Printing function for debugging.
    friend std::ostream& operator <<(std::ostream& os, const RivetONNXrt& rort){
      os << "RivetONNXrt Network Summary: \n";
      for (size_t i=0; i < rort._inNames.size(); ++i) {
        os << "- Input node " << i << " name: " << rort._inNames[i];
        os << ", dimensions: (";
        for (size_t j=0; j < rort._inDims[i].size(); ++j){
          if (j)  os << ", ";
          os << rort._inDims[i][j];
        }
        os << "), type (as ONNX enums): " << rort._inTypes[i] << "\n";
      }
      for (size_t i=0; i < rort._outNames.size(); ++i) {
        os << "- Output node " << i << " name: " << rort._outNames[i];
        os << ", dimensions: (";
        for (size_t j=0; j < rort._outDims[i].size(); ++j){
          if (j)  os << ", ";
          os << rort._outDims[i][j];
        }
        os << "), type (as ONNX enums): (" << rort._outTypes[i] << "\n";
      }
      return os;
    }

    /// Logger
    Log& getLog() const {
      string logname = "Rivet.RivetONNXrt";
      return Log::getLog(logname);
    }


  private:

    void getNetworkInfo() {

      Ort::AllocatorWithDefaultOptions allocator;

      // Retrieve network metadat
      _metadata = std::make_unique<Ort::ModelMetadata>(_session->GetModelMetadata());

      // find out how many input nodes the model expects
      const size_t num_input_nodes = _session->GetInputCount();
      _inDimsFlat.reserve(num_input_nodes);
      _inTypes.reserve(num_input_nodes);
      _inDims.reserve(num_input_nodes);
      _inNames.reserve(num_input_nodes);
      _inNamesPtr.reserve(num_input_nodes);
      for (size_t i = 0; i < num_input_nodes; ++i) {
        // retrieve input node name
        auto input_name = _session->GetInputNameAllocated(i, allocator);
        _inNames.push_back(input_name.get());
        _inNamesPtr.push_back(std::move(input_name));

        // retrieve input node type
        auto in_type_info = _session->GetInputTypeInfo(i);
        auto in_tensor_info = in_type_info.GetTensorTypeAndShapeInfo();
        _inTypes.push_back(in_tensor_info.GetElementType());
        _inDims.push_back(in_tensor_info.GetShape());
      }

      // Fix negative shape values - appears to be an artefact of batch size issues.
      for (auto& dims : _inDims) {
        int64_t n = 1;
        for (auto& dim : dims) {
          if (dim < 0)  dim = abs(dim);
          n *= dim;
        }
        _inDimsFlat.push_back(n);
      }

      // find out how many output nodes the model expects
      const size_t num_output_nodes = _session->GetOutputCount();
      _outDimsFlat.reserve(num_output_nodes);
      _outTypes.reserve(num_output_nodes);
      _outDims.reserve(num_output_nodes);
      _outNames.reserve(num_output_nodes);
      _outNamesPtr.reserve(num_output_nodes);
      for (size_t i = 0; i < num_output_nodes; ++i) {
        // retrieve output node name
        auto output_name = _session->GetOutputNameAllocated(i, allocator);
        _outNames.push_back(output_name.get());
        _outNamesPtr.push_back(std::move(output_name));

        // retrieve input node type
        auto out_type_info = _session->GetOutputTypeInfo(i);
        auto out_tensor_info = out_type_info.GetTensorTypeAndShapeInfo();
        _outTypes.push_back(out_tensor_info.GetElementType());
        _outDims.push_back(out_tensor_info.GetShape());
      }

      // Fix negative shape values - appears to be an artefact of batch size issues.
      for (auto& dims : _outDims) {
        int64_t n = 1;
        for (auto& dim : dims) {
          if (dim < 0)  dim = abs(dim);
          n *= dim;
        }
        _outDimsFlat.push_back(n);
      }
    }

  private:

    /// ONNXrt environment for this session
    std::unique_ptr<Ort::Env> _env;

    /// ONNXrt session holiding the network
    std::unique_ptr<Ort::Session> _session;

    /// Network metadata
    std::unique_ptr<Ort::ModelMetadata> _metadata;

    /// Input/output node dimensions
    ///
    /// @note Each could be a multidimensional tensor
    vector<vector<int64_t>> _inDims, _outDims;

    /// Equivalent length for flattened input/ouput node structure
    vector<int64_t> _inDimsFlat, _outDimsFlat;

    /// Types of input/output nodes (as ONNX enums)
    vector<ONNXTensorElementDataType> _inTypes, _outTypes;

    /// Pointers to the ONNXrt inout/output node names
    vector<Ort::AllocatedStringPtr> _inNamesPtr, _outNamesPtr;

    /// C-style arrays of the input/output node names
    vector<const char*> _inNames, _outNames;
  };

  /// Useful function for getting onnx file paths
  /// Based on getDatafilePath from RivetYODA.cc
  string getONNXFilePath(const string& filename) {
    /// Try to find an ONNX file matching this analysis name
    const string path1 = findAnalysisRefFile(filename);
    if (!path1.empty()) return path1;
    throw Rivet::Error("Couldn't find a ref data file for '" + filename +
                       "' in data path, '" + getRivetDataPath() + "', or '.'");
  }

  /// Function to get a RivetONNXrt object from an analysis name
  /// Use suffix to help disambiguate if an analysis requires 
  /// multiple networks.
  /// @todo: If ONNX is ever fully integrated into rivet, move
  /// to analysis class.
  unique_ptr<RivetONNXrt> getONNX(const string& analysisname, const string& suffix = ".onnx"){
    return make_unique<RivetONNXrt>(getONNXFilePath(analysisname+suffix));
  }
}
#endif
