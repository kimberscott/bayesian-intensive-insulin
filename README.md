# Setup instructions using pip
1. Clone the repo
2. Create and activate a virtual environment:
   A. Using virtualenv:
      * Install virtualenv if needed

        `pip3 install virtualenv`
      * Create and activate a virtual environment for this project

        `bayesian_intensive_insulin$ virtualenv -p python3 venv`
        `bayesian_intensive_insulin$ source venv/bin/activate`
        
   B. Using pyenv-virtualenv:
  
      * Install readline, xz, pyenv, and pyenv-virtualenv. Instructions below for OSX; see 
        [Pyenv docs](https://akrabat.com/creating-virtual-environments-with-pyenv/) for other OS.
            
        `brew install readline xz`
        `brew install pyenv pyenv-virtualenv` 
        
      * Install Python 3.10.8 using pyenv:
        
        `pyenv install 3.10.8`
    
      * Create a virtual environment using this version: 
    
        `pyenv virtualenv 3.10.8 insulin`
    
      * Automatically activate this environment when entering this directory: 
    
        `bayesian-intensive-insulin$ echo 'insulin' > .python-version`
  
3. From within your virtual environment, install all requirements:

   `bayesian_intensive_insulin$ pip install -r requirements.txt`
   `bayesian_intensive_insulin$ pip install -e .`
    
4. Add your environment as a kernel accessible to Jupyter:

   `ipython kernel install --user --name insulin`

4. Run `jupyter notebook`. It should open your browser and let you select `pymongo_diabetes.ipynb`
5. As you work, install any additional packages you need in your virtual environment; if you do, also 
   update the `requirements.txt` file.

# Making changes more readable: stripping down the .ipynb metadata

We would like to share changes to the code, but don't necessarily want to keep track of 
changes to the number of times each block has been executed, the output (actual blood sugar
values for particular weeks), etc. A good overview of the problem and a proposed solution 
using JQ is here: http://timstaley.co.uk/posts/making-git-and-jupyter-notebooks-play-nice/
To follow that approach:

* Install JQ following these directions: https://stedolan.github.io/jq/download/

* Test it out: the following command should output a stripped-down version of your notebook

    q --indent 1 \
        '
        (.cells[] | select(has("outputs")) | .outputs) = []
        | (.cells[] | select(has("execution_count")) | .execution_count) = null
        | .metadata = {"language_info": {"name":"python", "pygments_lexer": "ipython3"}}
        | .cells[].metadata = {}
        ' pymongo_diabetes.ipynb
    
* Follow the directions at http://timstaley.co.uk/posts/making-git-and-jupyter-notebooks-play-nice/ for 
editing your `.gitconfig` and `.gitattributes` files so that this is automatically run to 
clean your notebook upon committing changes.

* If you are using a git GUI interface like Sourcetree, note that it may define its own PATH and 
not be able to access `jq`. In this case, substitute the actual path to `jq` (e.g., `/usr/local/bin/jq`)
in the `clean` command in your `.gitconfig` file.

# TODO:

Setting up appropriate model assumptions:
* Allow varying breakfast/lunch/dinner carb ratios and sensitivity factors
* Look into using lognormal distribution for sigmaBG - does this matter for assessing certainty?
* Use more realistic priors for basal, CR, CF
* Use actual Humalog & carb action time profiles instead of linear
* Allow action time to scale

Accounting for known noise / deviations from model assumptions:
* Treat sugar to treat differently (expect to act faster, and/or just discard glucose measurements in following N hours)
* Treat exercise snacks differently (expect to have smaller influence, and/or discard glucose measurements in following N hours)
* Expect appropriate noise on insulin & carbs, and infer "true" IOB/COB profiles
* Allow complete discounting of a timepoint with some small probability - the "hrmmmm" factor for typos, REALLY off measurements, etc.

Basal profiles & pump data:
* Infer basal profile (in N-hour blocks) rather than overall rate
* Import pump basal data (incl. temp rates, basal IQ)

Evaluation and sharing:
* Evaluate model fit to data
* Conduct comparisons to actual protocols
* Code cleanup for readability
* Code cleanup for efficiency
* Set up to import other datasets