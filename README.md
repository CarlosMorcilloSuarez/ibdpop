ibdpop
                a tool to analyse IBD segments derived from individuals
                belonging to different populations.
        
        SYNOPSIS
        
            ibdpop <command> [options]
            
        COMMANDS
        
            version            Shows the current version of the program 
            
            help               Shows this help page
        
            matrixcounts       Matrix of the counts of IBDs between all pairs
                               of individuals
                               
            quality            Plots scores vs length of IBDs for QC
            
            recode             Creates a new dataset applying the current 
                               filters
            
        
        OPTIONS
                    
            --file <prefix>
                    Prefix of the files with IBDs (prefix.ibd) and individuals 
                    (prefix.ind) information
                    
            --out <prefix>
                    Prefix that will be used in all files created during the
                    execution
                    
            --color-code <file>
                    File containing the color codes for each population to be
                    used in plots.
                    
            --min-length <length>
                    Filters out IBD segments shorter of <length> cM
                    
            --impute-length
                    Imputes IBD segments length from chromosomal position 
                    instead of using provided cM length
                    
            --min-score <score>
                    Filters out IBD segments with score smaller that <score>
                    
            --dpi <dpi>
                    dpi (dots per inch) to be used plots
                    Default = 100
            
