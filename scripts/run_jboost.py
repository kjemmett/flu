import time
import os
import sys

# set classpath just in case
JBOOST_DIR = '/Users/kje/work/proj/flu/src/jboost'
os.environ['JBOOST_DIR'] = JBOOST_DIR 
try:
    os.environ['CLASSPATH'] = os.environ['CLASSPATH'] + ':' + ':'.join([JBOOST_DIR+'/dist/jboost.jar', JBOOST_DIR+'/lib/jfreechart-1.0.10.jar', JBOOST_DIR+'/lib/jcommon-1.0.8.jar']) 
except:
    os.environ['CLASSPATH'] = ':'.join([JBOOST_DIR+'/dist/jboost.jar', JBOOST_DIR+'/lib/jfreechart-1.0.10.jar', JBOOST_DIR+'/lib/jcommon-1.0.8.jar']) 

segment = 'concat.PB2.PB1.PA.HA.NP.NA.M1.NS1'

booster = 'AdaBoost'
rounds = 10
atreeoption = 'ADD_ALL'

jobname = '.'.join([segment, booster, str(rounds), atreeoption])
jobfolder = '.'.join([jobname, time.strftime("%m-%d-%y-%H-%M-%S")])

projpath = '/Users/kje/work/proj/flu/data/jboost'
jobpath = os.path.join(projpath, jobfolder)
specfile = os.path.join(projpath, segment+'.spec')
trainfile = os.path.join(projpath, segment+'.full')
testfile = os.path.join(projpath, segment+'.h7n9.test')

logfile = os.path.join(jobpath, jobname+'.log')

command = 'java -Xmx1000M -cp ' + os.getenv('CLASSPATH') \
        + ' jboost.controller.Controller -b ' + booster \
        + ' -p 3' \
        + ' -ATreeType ' + atreeoption \
        + ' -numRounds ' + str(rounds) \
        + ' -S ' + jobname \
        + ' -n ' + specfile \
        + ' -t ' + trainfile\
        + ' -T ' + testfile \
        + ' -a -1' \
        + ' -log ' + logfile \

# move to run folder
try:
    os.makedirs(jobpath)
except:
    pass

os.chdir(jobpath)

# run learner
error = os.system(command)
if (error != 0):
    sys.exit(1)


# plot output
plot_script = os.path.join(JBOOST_DIR,'scripts/atree2dot2ps.pl')
treefile = jobname+'.output.tree'
infofile = jobname+'.info'
plot_command = plot_script \
        + ' -i ' + infofile \
        + ' -s ' + '../'+os.path.basename(specfile) \
        + ' -t ' + treefile

error = os.system(plot_command)
if (error != 0):
    sys.exit(1)
