{gsub(/---[-]+/,"")}
{l[NR]=$$0}                   
END {                        
  for (i=1; i<=NR; i++)     
    if (l[i+1] ~ /^  "[^"]*"$/) {   
      d=l[i+1]                      
      sub(/^  "/,"",d); sub(/"$/,"",d)
      print "\n# " d "\n" l[i++]        
    } else 
      print l[i]  }               
