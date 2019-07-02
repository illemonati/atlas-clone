#include "TTaskList.h"

void TTaskList::fillTaskMap_all(){
    //copy regular tasks
    for(std::map< std::string, TTask* >::iterator it = taskMap_regular.begin(); it != taskMap_regular.end(); ++it){
        taskMap_all.insert(std::pair<std::string, TTask*>(it->first, it->second));
    }

    //copy debug tasks
    for(std::map< std::string, TTask* >::iterator it = taskMap_debug.begin(); it != taskMap_debug.end(); ++it){
        taskMap_all.insert(std::pair<std::string, TTask*>(it->first, it->second));
    }
}

TTask *TTaskList::getTask(const std::string task, TLog *logfile){
    std::map< std::string, TTask* >::iterator it = taskMap_all.find(task);
    if(it == taskMap_all.end()){
        logfile->setVerbose(true);
        printAvailableTasks(logfile);
        throwErrorUnknownTask(task);
    }
    return it->second;
}

void TTaskList::throwErrorUnknownTask(const std::string &Task){
    //calculate Levenshtein Distance and return closes match
    //initialize crazy values
    std::string bestMatch;
    unsigned int minDistance = 99999999;

    //do comparison iresspective of case: do all in lower case
    std::string task = Task;
    std::transform(task.begin(), task.end(), task.begin(), ::tolower);

    //now loop over all others
    for(std::map< std::string, TTask* >::iterator it = taskMap_regular.begin(); it != taskMap_regular.end(); ++it){
        //get string in lower case
        std::string s = it->first;
        std::transform(s.begin(), s.end(), s.begin(), ::tolower);

        //now compare
        unsigned int dist = levenshteinDistance(task, s, 0, 1); //insertion, deletion costs = 1,2
        if(dist < minDistance){
            bestMatch = it->first;
            minDistance = dist;
        }
    }

    //report best unless it has less than two characters overlapping
    if(minDistance < 2*(task.length()-1))
        throw "Unknown task '" + Task + "'! Did you mean '" + bestMatch  + "'?";
    else
        throw "Unknown task '" + Task + "'!";
}

TTaskList::TTaskList(){
    fillTaskMaps(taskMap_regular, taskMap_debug);

    //fill taskMap_all
    fillTaskMap_all();
}

TTaskList::~TTaskList(){
    for(std::map< std::string, TTask* >::iterator it = taskMap_all.begin(); it != taskMap_all.end(); ++it){
        delete it->second;
    }
}

bool TTaskList::taskExists(const std::string task){
    if(taskMap_all.find(task) != taskMap_all.end())
        return true;
    else
        return false;
}

void TTaskList::run(TParameters &parameters, TLog *logfile){
    std::string task = parameters.getParameterString("task");
    run(task, parameters, logfile);
}

void TTaskList::run(const std::string taskName, TParameters &parameters, TLog *logfile){
    TTask* task = getTask(taskName, logfile);
    task->run(taskName, parameters, logfile);
}

void TTaskList::printAvailableTasks(TLog *logfile){
    //get longest name
    unsigned int maxLen = 0;
    for(std::map< std::string, TTask* >::iterator it = taskMap_regular.begin(); it != taskMap_regular.end(); ++it){
        if(maxLen < it->first.length())
            maxLen = it->first.length();
    }

    //now print all tasks
    logfile->startNumbering("Available tasks:");
    for(std::map< std::string, TTask* >::iterator it = taskMap_regular.begin(); it != taskMap_regular.end(); ++it){
        std::string s = it->first;
        s.append(maxLen - it->first.length() + 3, ' '); //add white spaces to align explanations
        logfile->number(s + it->second->explanation());
    }

    logfile->endIndent();
}
