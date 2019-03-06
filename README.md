Bayesian causal modeling of diabetes data to infer insulin sensitivity and carb factors

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
