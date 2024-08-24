# Dacon competition
## About competition
IRAK4 IC50 활성 예측 모델 개발

## URL
https://dacon.io/competitions/official/236336/overview/description

## Notice
"_" 가 붙은 스크립트들은 데이터 viewing, 간단한 작업 위주의 내용입니다.   
대부분의 기능들은 .py 로 구현 예정입니다.   
현 시점에서는 rdkit 으로 SMILES 로부터 데이터 추출 완료 상태이며, 상세 내용은 Strategy 이하 내용 참조 바랍니다.  

# Strategy
After step 1, iterate step 2 and 3.

## Step 1. Extract features from SMILES code
※ Data extracted from SMILES are in "data/preprocessed" with their preprocessing method naming.  
Candidiates  
### [1] RDkit    
- Main page -> https://www.rdkit.org/
- Github -> https://github.com/rdkit/rdkit
- Document -> https://www.rdkit.org/docs/GettingStartedInPython.html  

#### 1. Morgan Fingerprints (Circular Fingerprints)   
- 설명:  
  - Morgan fingerprints는 화합물의 구조를 기반으로 한 분자 지문  
  - 특정 원자에서 출발하여 주변의 화학적 환경을 반경별로 해시 코드로 변환   
  - 반경이 커질수록 더 많은 화학적 특성이 포함  
- 장점: 고유한 화학적 환경을 효과적으로 캡처할 수 있어 분자 유사성 검색에 적합
- 도구: RDKit 
- DOC: https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-fingerprints

#### 2. MACCS Keys  
- 설명: MACCS (Molecular ACCess System) keys는 166개의 이진 템플릿을 사용하여 분자의 존재 여부를 표시하는 고정 길이의 이진 벡터를 생성  
- 장점: 특정 화학적 구조 모티프의 유무를 쉽게 캡처 가능. (간단하고 빠른 계산)  
- 도구: RDKit  
- DOC: https://www.rdkit.org/docs/GettingStartedInPython.html#maccs-keys

#### 3. Molecular Descriptors
- 설명: 분자 지표(Molecular Descriptors)는 분자의 다양한 물리적, 화학적 특성을 수치적으로 나타냄. 예를 들어, 분자량, 극성 표면적, 수소 결합 수 등이 존재.  
- 장점: 분자의 전반적인 물리적, 화학적 특성을 파악하는 데 유용하며, 머신러닝 모델의 성능을 향상 가능  
- 도구: RDKit 및 ChemAxon의 Marvin, Mordred 등의 툴에서 다양한 분자 지표를 계산  (RDKit 우선 사용)
- DOC: https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-descriptors

### [2]   
### [3]   



## Step 2. Select features



## Step 3. Modeling & Evaluation 
[1] Validation: 10-fold  
Train:Validation:Test = 7:2:1  

[2] Machine-Learning based (maybe tree-based boosting?)  
Basic: Xgboost, Lgboost, CatBoost ...   
What else: ?

[3] Neural-network architecture  
Custom:   
Pre-built: 

  
  

# Settings
## How to set virtual environment in your local machine
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

