# scBC

There are the code and data of A protocol of using scBC to analysis FGM perturbation during the progression of diseases, and a more detailed description of the properties and methods of scBC can be found here.

## Installation

To install our package, run

```bash
conda install -a conda-forge scvi-tools #pip install scvi-tools
pip install scBC
```

## Quick start

To simply illustrate how one can run our scBC model, here we use the subsampled HEART dataset as an example. 

```python
from scBC import data
data = ad.read_h5ad("[your_h5ad_file].h5ad")
```

If you can't download the data automatically, you can also download them manually from our [repository](https://github.com/GYQ-form/scBC/tree/main/data) and place them in your working directory.

Now you are ready to set up the scBC model:

```python
from scBC.model import scBC
my_model = scBC(adata=data,batch_key='cell_source')
```

Notice that `adata` must be an [AnnData](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html#anndata.AnnData) object. `batch_key` specify the specific column in adata.var to be used as batch annotation. 

We first train the VAE model with the default setting:

```python
my_model.train_VI()
```

Then we can get the reconstructed data, we use 10 posterior samples to estimate the mean expression:

```python
my_model.get_reconst_data(n_samples=10)
```

Since scBC encourage prior information about genes' co-expression to be provided,  before conducting biclustering, we first generate a prior edge. Thankfully, we have also provide an API for automatically generating prior edge from biomart database, just run with default setting:

```python
my_model.get_edge()
```

Now a edge array has been stored in my_model.edge. Run biclustering will automatically use it, we set the number of biclusters as 3 here:

```python
my_model.Biclustering(L=3)
```

Take a look at the results:

```python
my_model.S
```

Or maybe we are more interest in the strong classification result at cell level:

```python
my_model.strong_class()
```



## Parameter lists

Here we give some parameter lists for readers to use flexibly:

### model

#### scBC.model.scBC()

- **adata**:    *AnnData object with an observed n $\times$ p matrix, p is the number of measurements and n is the number of subjects.*  
- **layer**:    *if not None, uses this as the key in adata.layers for raw count data.*
- **batch_key**:     *key in adata.obs for batch information. Categories will automatically be converted into integer categories and saved to adata.obs['_scvi_batch']. If None, assigns the same batch to all the data.* 
- **labels_key**:    *key in adata.obs for label information. Categories will automatically be converted into integer categories and saved to adata.obs['_scvi_labels']. If None, assigns the same label to all the data.*
- **size_factor_key**:    *key in adata.obs for size factor information. Instead of using library size as a size factor, the provided size factor column will be used as offset in the mean of the likelihood. Assumed to be on linear scale.*
- **categorical_covariate_keys**:    *keys in adata.obs that correspond to categorical data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects  on the latent space). Thus, these should not be used for biologically-relevant factors that you do not want to correct for.*
- **continuous_covariate_keys**:  *keys in adata.obs that correspond to continuous data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space).  Thus, these should not be used for biologically-relevant factors that you do not want to correct for.*

#### Attributes

- `self.adata` - The AnnData used when initializing the scBC object.
- `self.p` - Number of  measurements(genes) of the expression matrix
- `self.n` - Number of  subjects(cells) of the expression matrx
- `self.vi_model` - The scvi model
- `self.reconst_data` - Reconstructed expression matrix with n $\times$ p
- `self.W` - **W** matrix after biclustering
- `self.Z` - **Z** matrix after biclustering
- `self.mu` - Parameter matrix $\pmb\mu$ after biclustering
- `self.convg_trail` - The converge trail of biclustering process (likelihood value)
- `self.S` - A list containing L dictionaries. The $i_{th}$ dictionary is the $i_{th}$ bicluster's results  
- `self.niter` - Iteration time during biclustering
- `self.edge` - The prior edge array

---

### Methods

#### scBC.train_VI()

- **n_hidden**:    Number of nodes per hidden layer.
- **n_latent**:    Dimensionality of the latent space.
- **n_layers**:    Number of hidden layers used for encoder and decoder NNs.
- **dropout_rate**:    Dropout rate for neural networks.
- **dispersion**:    One of the following:
  -  ``'gene'`` - dispersion parameter of NB is constant per gene across cells
  -  ``'gene-batch'`` - dispersion can differ between different batches
  -  ``'gene-label'`` - dispersion can differ between different labels
  -  ``'gene-cell'`` - dispersion can differ for every gene in every cell

- **gene_likelihood**:     One of:
  -  ``'nb'`` - Negative binomial distribution
  -  ``'zinb'`` - Zero-inflated negative binomial distribution
  -  ``'poisson'`` - Poisson distribution

- **latent_distribution**:    One of:
  -  ``'normal'`` - Normal distribution
  -  ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)

- **max_epochs**:    Number of passes through the dataset. If `None`, defaults to `np.min([round((20000 / n_cells) * 400), 400])`
- **use_gpu**:    Use default GPU if available (if None or True), or index of GPU to use (if int), or name of GPU (if `str`, e.g., `'cuda:0'`), or use CPU (if False).
- **batch_size**:    Minibatch size to use during training.
- **early_stopping**:    Perform early stopping.

---

#### scBC.get_reconst_data()

- **n_samples**:    Number of posterior samples to use for estimation.
- **batch_size**:    Minibatch size for data loading into model.

---

#### scBC.get_edge()

- **gene_list**:    A list containing HGNC gene names used to find prior edge in hsapiens_gene_ensembl dataset. Different actions will be performed based on the data type you provide:
  - `str` - take self.adata.var[gene_list] as the gene set
  - `list` - use the provided list as the gene set to be searched
  - `None(default)` --- use the variable names(self.adata.var_names) as the gene set

- **dataset**:    The dataset to be used for extracting prior information. Default is <u>hsapiens\_gene\_ensembl</u>.
- **intensity**:    A integer denoting the intensity of the prior. Generally speaking, the larger the number, the more prior information returned (the bigger the edge array is).

---

#### scBC.Biclustering()

- **L**:      number of the maximum biclusters
- **mat**:     a p x n numpy array used for biclustering. We will use the reconstructed data after calling scBC.get_reconst_data(). However, you can explicitly provide one rather than using the reconstructed  data.
- **edge**:     a 2-colounm matrix providing prior network information of some measurements. We prefer to use the edge given here, if not provided, we will use scBC.edge as a candidate. Edge can be automaticly generated using scBC.get_edge(). Running Biclustering() without any edge prior is also feasible.
- **dist_type**:  a p-element vector, each element represents the distribution type of measurements. 0 is Gaussian, 1 is Binomial, 2 is Negative Binomial, 3 is Poisson. If not given, all measurements are deemed to follow Gaussian distribution.
- **param**:    a p-element vector, each element is the required parameter for the correspondence distribution. $\zeta$ for Gaussian, $n_j$ for Binomial, $r_j$ for Negative Binomial, $N$ for Poisson. Default is set as a 1 vector. 
- **initWZ**:     select between "random" and "svd". If "random", randomly generated N(0,1) numbers are used, o.w. use SVD results for initialization.
- **maxiter**:   allowed maximum number of iterations.
- **tol**:     desired total difference of solutions.

---


