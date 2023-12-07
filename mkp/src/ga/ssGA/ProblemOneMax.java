/*
 * ProblemOneMax.java
 *
 * Created on 30 de octubre de 2001, 13:20
 */

package ga.ssGA;

/**
 *
 * @author  Antonio
 * @version 
 */
public class ProblemOneMax extends Problem{

    /** Creates new ProblemOneMax */
    //public ProblemOneMax() {
    //    super() ;
    //}

    public double Evaluate(Individual Indiv) {
      return ONEMAX(Indiv) ;
    }


  //    PRIVATE METHODS

  // Count the number of 1's in the string
  private double ONEMAX(Individual indiv)
  {
    double f=0.0;
    for(int i=0; i<CL; i++)
    if(indiv.get_allele(i)==1)
    f=f+1.0;
    indiv.set_fitness(f);
    return f;
  }

}
