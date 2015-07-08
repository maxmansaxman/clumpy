
import random

def randomNumberList(length, displayCounter = True):
    '''Generates a list of random numbers, at a length specified by the user.
    Additional option for whether a counter showing progress of generation
    (in percent) in displayed
    '''
    length = int(length)
    rng = random.SystemRandom()
    values = []
    counter = 0
    if displayCounter:
        # note: It's about 15% slower to display the progress than to not
        for i in range(length):
            values.append(rng.random())
            counter += 1
            if (counter * 100) % length == 0:
                print(str((counter*100)/length) + '% done')
    else:
        for i in range(length):
            values.append(rng.random())


    return values
