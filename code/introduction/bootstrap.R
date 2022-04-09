load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

library(dplyr)
library(xtable)

# das folgende ist so ähnlich wie im Skript, unter investigationg estimator quality by simulation
# ich kenne meine ursprüngliche Verteilung, hier Gleichverteilung auf 0,1 , Parameter bekannt, theta =1
# Statt des Parameters theta schätze ich seine Hälfte, und ich nehme das empirische Mittel
# Ich kann z.B. dafür benutzen:
# auf bias untersuchen
# Stichprobengröße bestimmen : wenn ich eine bestimmte sd für meinen Schätzer brauche, kann ich das n entsprechend wählen
# Schätzer vergleichen : oder wenn ich mehrere (unbiased) Schätzer habe kann ich sie bzgl ihrer sd vergleichen

# das sollte eigentlich genauso sein wie im Skript
f1 <- function(n){
  Vectorize(function(i){ 
    mean(runif(n))
  })(1:1000) %>% 
    (function(x) c(mean(x), sd(x)))
}
Vectorize(f1)(c(20,50,500))

# Änderung: anstatt ständig neue Daten zu generieren,
# werden nur einmal die zufälligen Daten generiert (a)
# dann nutze ich bootstraps (b)
f2 <- function(n){
  a <-  runif(n)
  Vectorize(function(i){ 
    b <- sample(a,length(a),T)
    mean(b)
  })(1:1000) %>% 
    (function(x) c(mean(x), sd(x)))
}
Vectorize(f2)(c(20,50,500))

# Jetzt muss ich das auf meine Situation übertragen:
# eigentlich habe ich nur eine Stichprobe
# meinen Cleveland Datensatz mit 303 observations (Patienten)

set.seed(100)

a <-  runif(303)
Vectorize(function(i){ 
  b <- sample(a,length(a),T)
  mean(b)
})(1:1000) -> res1
res1 %>% 
  (function(x) c(mean(x), sd(x)))

#res1 %>% density %>% plot
res1 %>% hist(main='this is what I can get from the Cleveland bootstraps\nfor a specified method'
              , xlab = ' logloss')

mean(a) # should be 0.5115 (and more generally, most likely not 0.5)

# Parallele zu Cleveand

# a = Cleveland
# b = bootstrapped Cleveland data with new CAD (the new CAD is why we call it simulation??)

# mean(b) is now a little more complicated:
# by applying a chosen method 
# create a small forest (using only Cleveland data. using OOB observations that were not sampled into b to rank trees by performance)
# and calculate its logloss on data.test,

# I should do this for the default forest, the regular small forest as base cases
# and for the new methods

# if this was a simulation as with the uniform random variables, I would compare to the true parameter.

# since I do not know the true parameter (true logloss)
# I could take the default forests as true parameter...
# maybe I should look at the logloss ratio (quotient) or the logloss overshoot (difference) for each method?
# that would look like this:

(res1 - mean(a)) %>% hist(main='comparing logloss for a specified method\nto logloss of default forest'
                          , xlab = 'logloss overshoot')
