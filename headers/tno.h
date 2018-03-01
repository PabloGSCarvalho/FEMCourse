/***************************************************************************
 *   Copyright (C) 2005 by Philippe R. B. Devloo                           *
 *   phil@corona                                                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef TNOS_H
#define TNOS_H


#include "TVec.h"

#include <fstream>
#include <iostream>

/**
Implementa um nó
@author Philippe R. B. Devloo
*/
class TNo{
public: // Construtores e destrutor
    /**
     * Construtor vazio - todas as variáveis são inicializadas com
     * valores nulos ou que indiquem falta de inicialização
     */
    TNo();
    
    /**
     * Construtor padrão - todos as variáveis são inicializadas 
     * com valores fornecidos
     */
    TNo (TVec<double> &coordenadas);
    
    /**
     * Construtor de cópia - as variáveis são inicializadas com valores 
     * do nó fornecido como parâmetro
     */
    TNo (const TNo &copy); 

    /**
     * Destrutor
     */
    ~TNo();
    
    
public: //Acesso a escrita dos objetos da classe
    /**
     * Define os valores dos dados da classe
     */
    void setData(TVec<double> &coordenadas);
    

    
public: //Acesso a leitura dos objetos da classe
    /**
     * Retorna a iésima coordenada
     */
    double &Co(int i);
       
    
  bool operator==(const TNo &other) const;
  
public: //Diversos
  /**
   * Imprime no local definido os dados da classe
   */
  void Print(std::ostream & out=std::cout) const;
  
  /**
   * Rotinas de teste
   */
   static int main();

protected:
    
    /**
     * Vetor para armazenar as coordenadas
     */
    double fCoordenadas[2];
    
};

inline std::ostream &operator<<(std::ostream &out, const TNo &obj)
{
  obj.Print(out);
  return out;
}

#endif
