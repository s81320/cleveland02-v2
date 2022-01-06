# follow
#https://github.com/llendway/github_for_collaboration/blob/master/github_for_collaboration.md

#install.packages('usethis')
library(usethis)
use_git_config(user.name = "silke"
               , user.email ="s81320@beuth-hochschule.de")

# set up a PAT (personal access token)
usethis::create_github_token()
# result : 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

#install.packages('gitcreds')
library(gitcreds)
gitcreds::gitcreds_set()
# if there are options given: replace current credentials and give the token above
# or first check the current password / token


