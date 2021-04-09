# SMYLE-analysis
A shared repository for jupyter notebooks examining SMYLE output

## Setting up
Clone this repo to your local directory, then create your own branch to work in:
```
git checkout -b <nameofyourbranch>
```

## Creating environment
```bash
conda env update -f envs/environment.yml
```

## Installing SMYLEutils (used in some of the example notebooks)
From within the local directory where you've cloned this repo and after activating your environment run
```bash
pip install -e . --user
```
Then when you're using this environment, all the modules in ./SMYLEutils will be available for you to import


