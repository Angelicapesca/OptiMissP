# OptiMissP Dashboard

The dashboard allows the qualitative and visual assessment of the missingness of data-independent acquisition mass spectrometry proteomic data. It contributes to define the better informed combination of imputation method and missingess threshold for the handling of missing values. 

The first panel allows the upload of the original not imputed proteomic dataset and its imputation. Otherwise, the user can upload the imputed dataset too. 

All the four panels are devoted to different types of data exploration:
* **Missingness' distribution**: two histrograms show the distribution of the number of missing values for instances (patients) and features (proteins)
* **Protein intensity distribution**: protein intensity distributions of imputed and not imputed data for different user selected missingness thresholds are put in comparison
* **Single protein intensity distribution**: chosen a specific protein, its protein distributions of imputed and not imputed data are comparatively plotted 
* **Topological Data Analysis of proteomic data**: Topological Data Analysis enables the creation of topologies from imputed and not imputed data for different missingness thresholds
