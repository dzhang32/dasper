##### Using sklearn in R via basilisk #####

env_sklearn <- basilisk::BasiliskEnvironment(
    envname = "env_sklearn",
    pkgname = "dasper",
    packages = c("scikit-learn==0.23.2", "pandas==0.25.3")
)
