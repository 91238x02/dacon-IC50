# Dacon competition
## About competition
IRAK4 IC50 활성 예측 모델 개발

## URL
https://dacon.io/competitions/official/236336/overview/description

  
  

# Strategy
After step 1, iterate step 2 and 3.

## Step 1. Extract features from SMIELS code
※ Data extracted from SMILES are in "data/preprocessed" with their preprocessing method naming.  
Candidiates  
[1] RDkit    
[2]   
[3]   



## Step 2. Select features



## Step 3. Training models with "train.csv" & evaluation 
[1] Validation: 10-fold  
Train:Validation:Test = 7:2:1  

[2] Machine-Learning based (maybe tree-based boosting?)  
Basic: Xgboost, Lgboost, CatBoost ...   
What else: ?

[3] Neural-network architecture  
Custom:   
Pre-built: 

  
  

# Settings
## How to set virtual environmet in your local machine
※ Python version: 3.11.8  
※ Commands below are for windows user  
$ cd {YOUR_PROJECT_DIR_PATH}  
$ git clone https://github.com/91238x02/dacon-IC50.git  
$ python -m venv .venv   
$ .venv\Scripts\activate   
$ pip install -r requirements.txt  

## DL Framework 
Pytorch -> v2.3.1  
$ pip install torch==2.3.1 torchvision==0.18.1 torchaudio==2.3.1 --index-url https://download.pytorch.org/whl/cu118  

## GPU
CUDA 11.8  
CUDNN 8.7  

