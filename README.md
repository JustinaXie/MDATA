## M-DATA: Multi-trait for De novo Mutation Association Test with Annotations
*Yuhan Xie, Mo Li, Weilai Dong, Wei Jiang, Hongyu Zhao*

To install the software for *M-DATA*, you need package *devtools* first

```{R}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
library(devtools)
```

Then, you can install *M-DATA* 
```{R}
install_github("JustinaXie/MDATA")
```

For software tutorial, please see https://github.com/JustinaXie/MDATA/blob/main/Vignette-%20User's%20Guide.Rmd


For detailed introduction, please refer to https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009849

## Credits
Those who use M-DATA should cite [M-DATA: A Statistical Approach to Jointly Analyzing De Novo Mutations for Multiple Traits](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009849)
