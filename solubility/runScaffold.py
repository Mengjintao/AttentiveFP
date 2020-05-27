import os
import torch
import torch.autograd as autograd
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.utils.data as Data
from scaffold import scaffold_split 

torch.manual_seed(8)

import time
import numpy as np
import gc
import sys
sys.setrecursionlimit(50000)
import pickle
torch.backends.cudnn.benchmark = True
torch.set_default_tensor_type('torch.cuda.FloatTensor')
# from tensorboardX import SummaryWriter
torch.nn.Module.dump_patches = True
import copy
import pandas as pd
#then import my own modules
from AttentiveFP import Fingerprint, Fingerprint_viz, save_smiles_dicts, get_smiles_dicts, get_smiles_array, moltosvg_highlight


from rdkit import Chem
# from rdkit.Chem import AllChem
from rdkit.Chem import QED
from rdkit.Chem import rdMolDescriptors, MolSurf
from rdkit.Chem.Draw import SimilarityMaps
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
#%matplotlib inline
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.cm as cm
import matplotlib
import seaborn as sns; sns.set_style("darkgrid")
from IPython.display import SVG, display
import sascorer
import itertools
from sklearn.metrics import r2_score
import scipy
import sys
from sklearn.utils import shuffle


task_name = 'solubility'
tasks = ['logS', 'weight']

def setup(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True

def train(model, dataset, optimizer, loss_function, batch_size, smiles_tasks_df, feature_dicts, epoch):
    model.train()
    np.random.seed(epoch)
    valList = np.arange(0,dataset.shape[0])
     
    #shuffle them
    np.random.shuffle(valList)
    batch_list = []
    for i in range(0, dataset.shape[0], batch_size):
        batch = valList[i:i+batch_size]
        batch_list.append(batch)   
    for counter, train_batch in enumerate(batch_list):
        batch_df = dataset.loc[train_batch,:]
        smiles_list = batch_df.cano_smiles.values
        y_val = batch_df[tasks[0]].values
        if (tasks[1] in smiles_tasks_df.columns.values):
            y_wgt = batch_df[tasks[1]].values
        else:
            y_wgt = torch.ones(y_val.shape)

        x_atom, x_bonds, x_atom_index, x_bond_index, x_mask, smiles_to_rdkit_list = get_smiles_array(smiles_list,feature_dicts)
        atoms_prediction, mol_prediction = model(torch.Tensor(x_atom),torch.Tensor(x_bonds),torch.cuda.LongTensor(x_atom_index),torch.cuda.LongTensor(x_bond_index),torch.Tensor(x_mask))
        
        model.zero_grad()
        loss = loss_function(mol_prediction, torch.Tensor(y_val).view(-1,1))    
        loss = loss * torch.Tensor(y_wgt).view(-1,1) 
        loss = loss.mean()
        loss.backward()
        optimizer.step()

def eval(model, dataset, batch_size, smiles_tasks_df, feature_dicts):
    model.eval()
    test_MAE_list = []
    test_MSE_list = []
    valList = np.arange(0,dataset.shape[0])
    batch_list = []
    for i in range(0, dataset.shape[0], batch_size):
        batch = valList[i:i+batch_size]
        batch_list.append(batch) 
    for counter, test_batch in enumerate(batch_list):
        
        batch_df = dataset.loc[test_batch,:]
        smiles_list = batch_df.cano_smiles.values
        y_val = batch_df[tasks[0]].values
        if (tasks[1] in smiles_tasks_df.columns.values):
            y_wgt = batch_df[tasks[1]].values
        else:
            y_wgt = torch.ones(y_val.shape)

        x_atom, x_bonds, x_atom_index, x_bond_index, x_mask, smiles_to_rdkit_list = get_smiles_array(smiles_list,feature_dicts)
        atoms_prediction, mol_prediction = model(torch.Tensor(x_atom),torch.Tensor(x_bonds),torch.cuda.LongTensor(x_atom_index),torch.cuda.LongTensor(x_bond_index),torch.Tensor(x_mask))
        MAE = F.l1_loss(mol_prediction, torch.Tensor(y_val).view(-1,1), reduction='none')        
        MSE = F.mse_loss(mol_prediction, torch.Tensor(y_val).view(-1,1), reduction='none')
        MAE = MAE * torch.Tensor(y_wgt).view(-1,1)
        MSE = MSE * torch.Tensor(y_wgt).view(-1,1)

#       print(x_mask[:2],atoms_prediction.shape, mol_prediction,MSE)
        se = MAE.data.squeeze().cpu().numpy()
        ae = MSE.data.squeeze().cpu().numpy()
          
        se = np.atleast_1d(se)
        ae = np.atleast_1d(ae)
#       test_MAE_list.extend(MAE.data.squeeze().cpu().numpy())
#       test_MSE_list.extend(MSE.data.squeeze().cpu().numpy())
        test_MAE_list.extend(ae)
        test_MSE_list.extend(se)
    return np.array(test_MAE_list).mean(), np.array(test_MSE_list).mean()

random_seed = 108 
batch_size = 128
epochs = 200
output_units_num = 1 # for regression model

def run(radius, T, fingerprint_dim, weight_decay, learning_rate, p_dropout, direction = False):
    radius = int (radius)
    T = int (T)
    fingerprint_dim = int (fingerprint_dim)
    print ("parameters")
    print (radius, T, fingerprint_dim, weight_decay, learning_rate, p_dropout)

    #raw_filename = "../data/delaney-processed.csv"
    #raw_filename = "~/jtmeng/SolCuration/org/esol/esol_org.csv"
    raw_filename = str(sys.argv[1])
    model_path = str(sys.argv[2])

    torch.cuda.set_device(int (sys.argv[3]))
    start_time = str(time.ctime()).replace(':','-').replace(' ','_')
    feature_filename = raw_filename.replace('.csv','.pickle')
    filename = raw_filename.replace('.csv','')
    prefix_filename = raw_filename.split('/')[-1].replace('.csv','')
    smiles_tasks_df = pd.read_csv(raw_filename)
    print (type(smiles_tasks_df))
    print (raw_filename)
    smilesList = smiles_tasks_df.smiles.values
    print("number of all smiles: ",len(smilesList))
    atom_num_dist = []
    remained_smiles = []
    canonical_smiles_list = []
    for smiles in smilesList:
        try:        
            mol = Chem.MolFromSmiles(smiles)
            atom_num_dist.append(len(mol.GetAtoms()))
            remained_smiles.append(smiles)
            canonical_smiles_list.append(Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True))
        except:
            print(smiles)
            pass
    print("number of successfully processed smiles: ", len(remained_smiles))
    smiles_tasks_df = smiles_tasks_df[smiles_tasks_df["smiles"].isin(remained_smiles)]
    # print(smiles_tasks_df)
    smiles_tasks_df['cano_smiles'] =canonical_smiles_list

    plt.figure(figsize=(5, 3))
    sns.set(font_scale=1.5)
    ax = sns.distplot(atom_num_dist, bins=28, kde=False)
    plt.tight_layout()
    # plt.savefig("atom_num_dist_"+prefix_filename+".png",dpi=200)
    plt.show()
    plt.close()

    if os.path.isfile(feature_filename):
        feature_dicts = pickle.load(open(feature_filename, "rb" ))
    else:
        feature_dicts = save_smiles_dicts(smilesList, filename)

    # feature_dicts = get_smiles_dicts(smilesList)
    remained_df  = smiles_tasks_df[smiles_tasks_df["cano_smiles"].isin(feature_dicts['smiles_to_atom_mask'].keys())]
    uncovered_df = smiles_tasks_df.drop(remained_df.index)
    print("not processed items")
    print (uncovered_df)

    seed = 5
    smiSet = remained_df["cano_smiles"].tolist()
    idxtrain, idxeval, idxtest = scaffold_split(smiSet, sizes=[0.8, 0.1, 0.1], balanced=True, seed=seed, logger=print)      
    print (idxtrain)
    print (idxeval)
    print (idxtest)

    all_scores = []
    for random_seed in range(5):
    #   remained_df.sample(n=len(remained_df), random_state=random_seed)    
        remained_df = remained_df.reset_index(drop=True)
        smiSet      = remained_df["cano_smiles"].tolist()
        
        idxtrain, idxeval, idxtest = scaffold_split(smiSet, sizes=[0.8, 0.1, 0.1], balanced=True, seed=random_seed, logger=print)
        train_df = remained_df.loc(idxtrain)
        valid_df = remained_df.loc(idxeval)
        test_df  = remained_df.loc(idxtest)

        valid_df = training_data.sample(frac=1/9, random_state=random_seed) # validation set
        train_df = training_data.drop(valid_df.index) # train set

        train_df = train_df.reset_index(drop=True)
        valid_df = valid_df.reset_index(drop=True)
        test_df  = test_df.reset_index(drop=True)

    # print(len(test_df),sorted(test_df.cano_smiles.values))

        x_atom, x_bonds, x_atom_index, x_bond_index, x_mask, smiles_to_rdkit_list = get_smiles_array([canonical_smiles_list[0]],feature_dicts)
        num_atom_features = x_atom.shape[-1]
        num_bond_features = x_bonds.shape[-1]
        loss_function = nn.MSELoss(reduction='none')
        model = Fingerprint(radius, T, num_atom_features, num_bond_features,fingerprint_dim, output_units_num, p_dropout)
        model.cuda()

        # optimizer = optim.Adam(model.parameters(), learning_rate, weight_decay=weight_decay)
        optimizer = optim.Adam(model.parameters(), 10**-learning_rate, weight_decay=10**-weight_decay)
        # optimizer = optim.SGD(model.parameters(), 10**-learning_rate, weight_decay=10**-weight_decay)

    # tensorboard = SummaryWriter(log_dir="runs/"+start_time+"_"+prefix_filename+"_"+str(fingerprint_dim)+"_"+str(p_dropout))

        model_parameters = filter(lambda p: p.requires_grad, model.parameters())
        params = sum([np.prod(p.size()) for p in model_parameters])
#        print(params)
#        for name, param in model.named_parameters():
#            if param.requires_grad:
#                print(name, param.data.shape)


        best_param ={}
        best_param["train_epoch"] = 0
        best_param["valid_epoch"] = 0
        best_param["train_MSE"] = 9e8
        best_param["valid_MSE"] = 9e8

        file_name = prefix_filename + '_' + str(radius) + '_' + str(T) + '_' + str(fingerprint_dim) + '_' + str(weight_decay) + '_' + str(learning_rate)
        for epoch in range(epochs):
            train_MAE, train_MSE = eval(model, train_df, batch_size, smiles_tasks_df, feature_dicts)
            valid_MAE, valid_MSE = eval(model, valid_df, batch_size, smiles_tasks_df, feature_dicts)
#     tensorboard.add_scalars('MAE',{'train_MAE':valid_MAE, 'test_MAE':valid_MSE}, epoch)
#     tensorboard.add_scalars('MSE',{'train_MSE':valid_MAE, 'test_MSE':valid_MSE}, epoch)
            if train_MSE < best_param["train_MSE"]:
                best_param["train_epoch"] = epoch
                best_param["train_MSE"] = train_MSE
            if valid_MSE < best_param["valid_MSE"]:
                best_param["valid_epoch"] = epoch
                best_param["valid_MSE"] = valid_MSE
#           if valid_MSE < 0.35:
                torch.save(model, model_path+'/model_'+file_name+'_'+start_time+'_'+str(epoch)+'.pt')
            if (epoch - best_param["train_epoch"] >8) and (epoch - best_param["valid_epoch"] >10):        
                break
            print(epoch, np.sqrt(train_MSE), np.sqrt(valid_MSE))
            train(model, train_df, optimizer, loss_function, batch_size, smiles_tasks_df, feature_dicts, epoch)

# evaluate model
        best_model = torch.load(model_path+'/model_'+file_name+'_'+start_time+'_'+str(best_param["valid_epoch"])+'.pt')     

        best_model_dict = best_model.state_dict()
        best_model_wts = copy.deepcopy(best_model_dict)

        model.load_state_dict(best_model_wts)
        (best_model.align[0].weight == model.align[0].weight).all()
        test_MAE, test_MSE = eval(model, test_df, batch_size, smiles_tasks_df, feature_dicts)
        print("best epoch:",best_param["valid_epoch"],"\n","test RMSE:",np.sqrt(test_MSE))
        all_scores.append(np.sqrt(test_MSE))

    all_scores = np.array(all_scores)
    print(f'fold cross validation')
    # Report scores for each fold
#    for fold_num, scores in enumerate(all_scores):
#        print(f'Seed {fold_num} ==> test rmse = {np.nanmean(scores):.6f}')

    mean_score, std_score = np.nanmean(all_scores), np.nanstd(all_scores)
    print(f'Overall test rmse = {mean_score:.6f} +/- {std_score:.6f}')
    return -1*mean_score 

radius = int(sys.argv[4])
T = int(sys.argv[5])
fingerprint_dim =  int(sys.argv[6])
weight_decay = float (sys.argv[7])
learning_rate = float (sys.argv[8])
run(radius=radius, T=T, fingerprint_dim=fingerprint_dim, weight_decay=weight_decay, learning_rate=learning_rate, p_dropout=0.2)
