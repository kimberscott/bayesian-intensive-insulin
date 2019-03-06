# Setup instructions using pip
* Clone the repo
* Install virtualenv if needed

    `pip3 install virtualenv`
* Create and activate a virtual environment for this project

    `bayesian_intensive_insulin$ virtualenv -p python3 venv`
    `bayesian_intensive_insulin$ source venv/bin/activate`
* Install requirements, including jupyter

    `bayesian_intensive_insulin$ pip3 install -r requirements.txt`
* Run `jupyter notebook`. It should open your browser and let you select `pymongo_diabetes.ipynb`
* As you work, install any additional packages you need from your virtual environment. If 
  you do so, also run `pip3 freeze > requirements.txt` and commit the `requirements.txt` file.

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

# TODO (can move to issues):

* Take start date/duration as params instead of segment, and also fetch range of protocol vals
* Plot inferred against in-use parameters

* Allow varying breakfast/lunch/dinner carb ratios and sensitivity factors
* Look into using lognormal distribution for sigmaBG - does this matter for assessing certainty?
* Use more realistic priors for basal, CR, CF
* Use actual Humalog & carb action time profiles instead of linear
* Allow action time to scale

* Expect appropriate noise on insulin & carbs, and infer "true" IOB/COB profiles
* Infer basal profile rather than overall rate

* Treat sugar to treat differently (expect to act faster, and/or just discard glucose measurements in following N hours)
* Treat exercise snacks differently (expect to have smaller influence, and/or discard glucose measurements in following N hours)

* Import pump basal data
* Allow complete discounting of a timepoint with some small probability - the "hrmmmm" factor for typos, REALLY off measurements, etc.

* Evaluate model fit to data
* Conduct comparisons to actual protocols
* Code cleanup for readability
* Code cleanup for efficiency
* Set up to import other datasets