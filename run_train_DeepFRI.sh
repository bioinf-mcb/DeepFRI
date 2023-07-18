#!/bin/bash
# Script for training DeepFRI model

# cuda libraries to run tf2 on gpu
# module load slurm gcc cuda/10.1.105_418.39 cudnn/v7.6.2-cuda-10.1

export CUDA_VISIBLE_DEVICES=0
export XLA_FLAGS=--xla_gpu_cuda_data_dir=/root/miniconda3/envs/mDF-mat/pkgs/cuda-toolkit

main_dir=/path/to/directory/with/training/data/
out_dir=/path/to/output/directory/

mkdir -p ${out_dir}

graph_conv_dims="512 512 512"
fully_connected_dims="1024"
graph_conv_layer=GraphConv
ontology_name='molecular_function'
ontology='mf'
cmap_thresh=10.0
data_dir=/${main_dir}/TFRecords/
cmap_data=UNIPROT
model_name=${out_dir}/DeepFRI-${cmap_data}_${graph_conv_layer}_gcd_$(echo $graph_conv_dims | tr ' ' '-')_fcd_${fully_connected_dims}_ca_${cmap_thresh}_${ontology}

annot_fn=${main_dir}/annotations/${ontology_name}.tsv


echo "Training ${ontology}..."

python train_DeepFRI.py \
    -gcd ${graph_conv_dims} \
    -fcd ${fully_connected_dims} \
    -l2 2e-5 \
    -lr 0.0002 \
    -gc ${graph_conv_layer} \
    -e 50  \
    -bs 64 \
    -ont ${ontology} \
    -lm trained_models/lstm_lm.hdf5 \
    --cmap_type ca \
    --cmap_thresh ${cmap_thresh} \
    --annot_fn ${annot_fn} \
    --model_name ${model_name} \
    --train_tfrecord_fn ${data_dir}/${ontology_name}/train/${cmap_data}_${ontology}_train \
    --valid_tfrecord_fn ${data_dir}/${ontology_name}/val/${cmap_data}_${ontology}_val
