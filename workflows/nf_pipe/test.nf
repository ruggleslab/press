#!/usr/bin/env nextflow run

process perlTask {
    """
    #!/usr/bin/perl

    print 'Hi there!' . '\n';
    """
}

process pythonTask {
    """
    #!/opt/homebrew/anaconda3/envs/main/bin/python

    x = 'Hello'
    y = 'world!'
    print("%s - %s" % (x,y))
    """
}

workflow {
    // perlTask()
    pythonTask()
}