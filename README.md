# plmvpaLite
Simplified matlab code for interfacing with liblinear and the Princeton MVPA toolbox and conducting PLR classification

Note: liblinear has two known quirks (12/03/23) - 1) performance is different between the windows and linux compilations; need to investigate. 2) performance is coded according to the first category encountered in the training set (e.g., between faces and scenes, if the first trial is a scene this becomes class A) - plmvpaLite tracks this information and re-codes so that the category labels and performance metrics are always in the order you specify them when setting up the prolem in the code
