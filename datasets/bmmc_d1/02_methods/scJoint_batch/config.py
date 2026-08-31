# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import torch
import os

class Config(object):
    def __init__(self):
        DB = '10x'
        self.use_cuda = True
        self.threads = 1

        if not self.use_cuda:
            self.device = torch.device('cpu')
        else:
            self.device = torch.device('cuda:0')
        
        if DB == '10x':
            # DB info
            self.number_of_class = 22
            self.input_size = 18353
            self.rna_paths = ['/path/to/data/BMMC_d1/scjoint2/rna1_scjoint.npz',\
                              '/path/to/data/BMMC_d1/scjoint2/rna2_scjoint.npz',\
                                '/path/to/data/BMMC_d1/scjoint2/rna3_scjoint.npz']
            self.rna_labels = ['/path/to/data/BMMC_d1/scjoint2/test1rna_celltype_scjoint.txt',\
                               '/path/to/data/BMMC_d1/scjoint2/test2rna_celltype_scjoint.txt',\
                                '/path/to/data/BMMC_d1/scjoint2/test3rna_celltype_scjoint.txt']		
            self.atac_paths = ['/path/to/data/BMMC_d1/scjoint2/atac1_gene_scjoint.npz',\
            					'/path/to/data/BMMC_d1/scjoint2/atac2_gene_scjoint.npz',\
                                    '/path/to/data/BMMC_d1/scjoint2/atac3_gene_scjoint.npz']
            self.atac_labels = [] #Optional. If atac_labels are provided, accuracy after knn would be provided.
            self.rna_protein_paths = []
            self.atac_protein_paths = []
            
            # Training config            
            self.batch_size = 256
            self.lr_stage1 = 0.01
            self.lr_stage3 = 0.01
            self.lr_decay_epoch = 20
            self.epochs_stage1 = 20
            self.epochs_stage3 = 20
            self.p = 0.8
            self.embedding_size = 64
            self.momentum = 0.9
            self.center_weight = 1
            self.with_crossentorpy = True
            self.seed = 1
            self.checkpoint = ''
        
        elif DB == "MOp":
            self.number_of_class = 21
            self.input_size = 18603
            self.rna_paths = ['data_MOp/YaoEtAl_RNA_snRNA_10X_v3_B_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_snRNA_10X_v3_A_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_snRNA_10X_v2_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_snRNA_SMARTer_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_scRNA_10X_v3_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_scRNA_10X_v2_exprs.npz',\
                                'data_MOp/YaoEtAl_RNA_scRNA_SMARTer_exprs.npz']
            self.rna_labels = ['data_MOp/YaoEtAl_RNA_snRNA_10X_v3_B_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_snRNA_10X_v3_A_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_snRNA_10X_v2_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_snRNA_SMARTer_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_scRNA_10X_v3_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_scRNA_10X_v2_cellTypes.txt',\
                                'data_MOp/YaoEtAl_RNA_scRNA_SMARTer_cellTypes.txt']
            self.atac_paths = ['data_MOp/YaoEtAl_ATAC_exprs.npz',\
                                'data_MOp/YaoEtAl_snmC_exprs.npz']
            self.atac_labels = ['data_MOp/YaoEtAl_ATAC_cellTypes.txt',\
                                'data_MOp/YaoEtAl_snmC_cellTypes.txt']
            self.rna_protein_paths = []
            self.atac_protein_paths = []
            
            # Training config            
            self.batch_size = 256
            self.lr_stage1 = 0.01
            self.lr_stage3 = 0.001
            self.lr_decay_epoch = 20
            self.epochs_stage1 = 10
            self.epochs_stage3 = 10
            self.p = 0.8
            self.embedding_size = 64
            self.momentum = 0.9
            self.center_weight = 20
            self.with_crossentorpy = True
            self.seed = 1
            self.checkpoint = '' 
            
        elif DB == "db4_control":
            self.number_of_class = 7 # Number of cell types in CITE-seq data
            self.input_size = 17668 # Number of common genes and proteins between CITE-seq data and ASAP-seq
            self.rna_paths = ['data/citeseq_control_rna.npz'] # RNA gene expression from CITE-seq data
            self.rna_labels = ['data/citeseq_control_cellTypes.txt'] # CITE-seq data cell type labels (coverted to numeric) 
            self.atac_paths = ['data/asapseq_control_atac.npz'] # ATAC gene activity matrix from ASAP-seq data
            self.atac_labels = ['data/asapseq_control_cellTypes.txt'] # ASAP-seq data cell type labels (coverted to numeric) 
            self.rna_protein_paths = ['data/citeseq_control_adt.npz'] # Protein expression from CITE-seq data
            self.atac_protein_paths = ['data/asapseq_control_adt.npz'] # Protein expression from ASAP-seq data
            
            # Training config            
            self.batch_size = 256
            self.lr_stage1 = 0.01
            self.lr_stage3 = 0.01
            self.lr_decay_epoch = 20
            self.epochs_stage1 = 20
            self.epochs_stage3 = 20
            self.p = 0.8
            self.embedding_size = 64
            self.momentum = 0.9
            self.center_weight = 1
            self.with_crossentorpy = True
            self.seed = 1
            self.checkpoint = '' 



            

