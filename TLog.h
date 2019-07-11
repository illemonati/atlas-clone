/*
 * TLog.h
 *
 *  Created on: Oct 17, 2010
 *      Author: wegmannd
 */

#ifndef TLOG_H_
#define TLOG_H_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sys/time.h>

class TLog{
private:
	bool isFile;
	std::string filename;
	bool isVerbose;
	int numIndent;
	std::string indent, indetOnlyTabs;
	bool printWarnings;
	long lastLineStartInFile;
	int numberingLevel;
	std::vector<int> numberingIndex;
	struct timeval start, end;

public:
	std::ofstream file;

    TLog();
    ~TLog();

    void close();
    void printTimeSinceStartOfProgam();
    void openFile(std::string Filename);
    void setVerbose(bool Verbose);
    bool verbose();
    void suppressWarings();
    void showWarings();
    void newLine();
    std::string getFilename();

	template<typename T>
	void add(T out){
		if(isFile) file << out;
		if(isVerbose) std::cout << out;
    }

	//---------------------------------------------------------
	//Indent
	//---------------------------------------------------------
    void fillIndentString();
    void addIndent(int n=1);
    void removeIndent(int n=1);
    void clearIndent();

	template<typename T>
	void startIndent(T out){
		list(out);
		addIndent();
    }

	template<typename T, typename U, typename V>
	void startIndent(T first, U middle, V last){
		list(first, middle, last);
		addIndent();
    }

	template<typename T, typename U, typename V, typename W, typename X>
	void startIndent(T first, U second, V third, W fourth, X fifth){
		list(first, second, third, fourth, fifth);
		addIndent();
    }

	template<typename T>
	void startIndentFlush(T out){
		listFlush(out);
		addIndent();
    }

	template<typename T>
	void endIndent(T out){
		list(out);
		removeIndent();
    }

	template<typename T, typename U, typename V>
	void endIndent(T first, U middle, V last){
		list(first, middle, last);
		removeIndent();
    }

    void endIndent();

	//---------------------------------------------------------
	//Numbering
	//---------------------------------------------------------
    void addNumberingLevel();
    void removeNumberingLevel();

	template<typename T>
	void startNumbering(T out){
		list(out);
		addNumberingLevel();
    }

    void endNumbering();

	template<typename T>
	void endNumbering(T out){
		number(out);
		removeNumberingLevel();
    }

	template<typename T>
	void number(T out){
		if(numberingLevel<0) addNumberingLevel();
		if(isFile){
			file << indetOnlyTabs << numberingIndex[numberingLevel] << ") ";
			if(numberingIndex[numberingLevel] < 10) file << ' ';
			file << out << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose){
			std::cout << indetOnlyTabs << numberingIndex[numberingLevel] << ") ";
			if(numberingIndex[numberingLevel] < 10) std::cout << ' ';
			std::cout << out << std::endl;
		}
		++numberingIndex[numberingLevel];
    }

	template<typename T>
	void numberWithIndent(T out){
		number(out);
		addIndent();
    }

	template<typename T>
	void numberFlush(T out){
		if(numberingLevel<0) addNumberingLevel();
		if(isFile){
			file << indetOnlyTabs << numberingIndex[numberingLevel] << ')';
			if(numberingIndex[numberingLevel] < 10) file << ' ';
			file << out << std::flush;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose){
			std::cout << indetOnlyTabs << numberingIndex[numberingLevel] << ')';
			if(numberingIndex[numberingLevel] < 10) std::cout << ' ';
			std::cout << out << std::flush;
		}
		++numberingIndex[numberingLevel];
    }

	template<typename T>
	void overNumber(T out){
		if(isFile){
			file.seekp(lastLineStartInFile);
			file << indetOnlyTabs << numberingIndex[numberingLevel] << ')';
			file << out << std::flush;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose){
			std::cout << '\xd' << indetOnlyTabs << numberingIndex[numberingLevel]-1 << ')';
			std::cout << out << std::flush;
			std::cout << out << std::flush;
		}
    }


	//---------------------------------------------------------
	//Write
	//---------------------------------------------------------
	template<typename T>
	void write(T out){
		if(isFile){
			file << out << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << out << std::endl;
    }

	template<typename T, typename V>
	void write(T first, V last){
		if(isFile){
			file << first << last << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << first << last << std::endl;
    }

	template<typename T, typename U, typename V>
	void write(T first, U middle, V last){
		if(isFile){
			file << first << middle << last << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << first << middle << last << std::endl;
    }
    void done();

	//---------------------------------------------------------
	//List
	//---------------------------------------------------------
	template<typename T>
	void list(T out){
		if(isFile){
			file << indent << out << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << indent << out << std::endl;
    }

	template<typename T, typename U>
	void list(T first, U last){
		if(isFile){
			file << indent << first << last << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << indent << first << last << std::endl;
    }

	template<typename T, typename U, typename V>
	void list(T first, U middle, V last){
		if(isFile){
			file << indent << first << middle << last << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << indent << first << middle << last << std::endl;
    }

	template<typename T, typename U, typename V, typename W, typename X>
	void list(T first, U second, V third, W fourth, X fifth){
		if(isFile){
			file << indent << first << second << third << fourth << fifth << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << indent << first << second << third << fourth << fifth << std::endl;
    }

	//---------------------------------------------------------
	//Conclude
	//---------------------------------------------------------
	template<typename T>
	void conclude(T out){
		std::string temp="";
		for(int i=0; i<=numIndent; ++i) temp+="   ";
		if(isFile){
			file << temp << "-> " << out << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << temp << "-> " << out << std::endl;
    }

	template<typename T, typename U>
	void conclude(T first, U last){
		std::string temp="";
		for(int i=0; i<=numIndent; ++i) temp+="   ";
		if(isFile){
			file << temp << "-> " << first << last << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << temp << "-> " << first << last << std::endl;
    }

	template<typename T, typename U, typename V>
	void conclude(T first, U middle, V last){
		std::string temp="";
		for(int i=0; i<=numIndent; ++i) temp+="   ";
		if(isFile){
			file << temp << "-> " << first << middle << last << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << temp << "-> " << first << middle << last << std::endl;
    }

	template<typename T, typename U, typename V, typename W>
	void conclude(T first, U middle, V secondlast, W last){
		std::string temp="";
		for(int i=0; i<=numIndent; ++i) temp+="   ";
		if(isFile){
			file << temp << "-> " << first << middle << secondlast <<  last << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << temp << "-> " << first << middle << secondlast << last << std::endl;
    }

	template<typename T, typename U, typename V, typename W, typename X>
	void conclude(T first, U second, V third, W fourth, X fifth){
		std::string temp="";
		for(int i=0; i<=numIndent; ++i) temp+="   ";
		if(isFile){
			file << temp << "-> " << first << second << third << fourth << fifth << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << temp << "-> " << first << second << third << fourth << fifth << std::endl;
    }

	//---------------------------------------------------------
	//overWrite
	//---------------------------------------------------------
	template<typename T>
	void overWrite(T out){
		if(isFile){
			file.seekp(lastLineStartInFile);
			file << out << std::endl;
			//file << '\xd'<< out << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << '\xd' << out << std::endl;
    }

	template<typename T, typename U, typename V>
	void overWrite(T first, U middle, V last){
		if(isFile){
			file.seekp(lastLineStartInFile);
			file << first << middle << last << std::endl;
			//file << '\xd'<< first << middle << last << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << '\xd' << first << middle << last << std::endl;
    }

	template<typename T, typename U, typename V, typename W, typename X>
	void overWrite(T first, U second, V third, W fourth, X fifth){
		if(isFile){
			file.seekp(lastLineStartInFile);
			file << first << second << third << fourth << fifth << std::endl;
			//file << '\xd'<< first << second << third << fourth << fifth << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << '\xd' << first << second << third << fourth << fifth << std::endl;
    }

	//---------------------------------------------------------
	//overList
	//---------------------------------------------------------
	template<typename T>
	void overList(T out){
		if(isFile){
			file.seekp(lastLineStartInFile);
			file << indent << out << std::endl;
			//file << '\xd'<< indent << out << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << '\xd' << indent << out << std::endl;
    }

	template<typename T, typename U, typename V>
	void overList(T first, U middle, V last){
		if(isFile){
			file.seekp(lastLineStartInFile);
			file << indent << first << middle << last << std::endl;
			//file << '\xd'<< indent << first << middle << last << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << '\xd' << indent << first << middle << last << std::endl;
    }

	template<typename T, typename U, typename V, typename W, typename X>
	void overList(T first, U second, V third, W fourth, X fifth){
		if(isFile){
			file.seekp(lastLineStartInFile);
			file << indent << first << second << third << fourth << fifth << std::endl;
			//file << '\xd'<< indent << first << second << third << fourth << fifth << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << '\xd' << indent << first << second << third << fourth << fifth << std::endl;
    }

	template<typename T>
	void writeFileOnly(T out){
		if(isFile){
			file << out << std::endl;
			lastLineStartInFile=file.tellp();
		}
    }

	template<typename T>
	void listFileOnly(T out){
		if(isFile){
			file << indent << out << std::endl;
			lastLineStartInFile=file.tellp();
		}
    }

	template<typename T>
	void listNoFile(T out){
		std::cout << indent << out << std::endl;
    }

	//---------------------------------------------------------
	//warning / error
	//---------------------------------------------------------
	template<typename T>
	void warning(T out){
		if(printWarnings){
			if(isFile){
				file << "WARNING: " << out << std::endl << std::endl;
				lastLineStartInFile=file.tellp();
			}
			std::cerr << "WARNING: " << out << std::endl << std::endl;
		}
    }

	template<typename T, typename U, typename V>
	void warning(T first, U middle, V last){
		if(printWarnings){
			if(isFile){
				file << "WARNING: " << first << middle << last << std::endl << std::endl;
				lastLineStartInFile=file.tellp();
			}
			std::cerr << "WARNING: " << first << middle << last << std::endl << std::endl;
		}
    }

	template<typename T>
	void error(T out){
		newLine();
		if(isFile){
			file << "ERROR: " << out << std::endl << std::endl;
			lastLineStartInFile=file.tellp();
		}
		std::cerr << "ERROR: " << out << std::endl << std::endl;
    }

	template<typename T, typename U, typename V>
	void error(T first, U middle, V last){
		newLine();
		if(isFile){
			file << "ERROR: " << first << middle << last << std::endl << std::endl;
			lastLineStartInFile=file.tellp();
		}
		std::cerr << "ERROR: " << first << middle << last << std::endl << std::endl;
    }

	//---------------------------------------------------------
	//Flush
	//---------------------------------------------------------
	template<typename T>
	void flush(T out){
			if(isFile) file << out << std::flush;
			if(isVerbose) std::cout << out << std::flush;
    }

	template<typename T, typename U, typename V>
	void flush(T first, U middle, V last){
		if(isFile) file << first << middle << last << std::flush;
		if(isVerbose) std::cout << first << middle << last << std::flush;
    }

	template<typename T>
	void listFlush(T out){
		if(isFile) file << indent << out << std::flush;
		if(isVerbose) std::cout << indent << out << std::flush;
    }

	template<typename T>
	void listFlush(T out, char symbol){
		if(isFile) file << indetOnlyTabs << symbol << ' ' << out << std::flush;
		if(isVerbose) std::cout << indetOnlyTabs << symbol << ' ' << out << std::flush;
    }

	template<typename T, typename U, typename V>
	void listFlush(T first, U middle, V last){
			if(isFile) file << indent << first << middle << last << std::flush;
			if(isVerbose) std::cout << indent << first << middle << last << std::flush;
    }

	template<typename T, typename U, typename V, typename W, typename X>
	void listFlush(T first, U second, V third, W fourth, X fifth){
			if(isFile) file << indent << first << second << third << fourth << fifth << std::flush;
			if(isVerbose) std::cout << indent << first << second << third << fourth << fifth << std::flush;
    }


	template<typename T>
	void flushFileOnly(T out){
		if(isFile) file << out << std::flush;
    }

	template<typename T>
	void listFlushFileOnly(T out){
		if(isFile) file << indent << out << std::flush;
    }

	template<typename T>
	void overFlush(T out){
			if(isFile){
				 file.seekp(lastLineStartInFile);
				 file << out << std::flush;
				 //file << '\xd' << out << std::flush;
			}
			if(isVerbose) std::cout << '\xd' << out << std::flush;
    }

	template<typename T, typename U, typename V>
	void overFlush(T first, U middle, V last){
		if(isFile){
			file.seekp(lastLineStartInFile);
			file << first << middle << last << std::flush;
			//file << '\xd' << first << middle << last << std::flush;
		}
		if(isVerbose) std::cout << '\xd' << first << middle << last << std::flush;
    }

	template<typename T, typename U, typename V, typename W, typename X>
	void overFlush(T first, U second, V third, W fourth, X fifth){
		if(isFile){
			file.seekp(lastLineStartInFile);
			file << first << second << third << fourth << fifth << std::flush;
			//file << '\xd' << first << second << third << fourth << fifth << std::flush;
		}
		if(isVerbose) std::cout << '\xd' << first << second << third << fourth << fifth << std::flush;
    }

	template<typename T>
	void listOverFlush(T out){
			if(isFile){
				file.seekp(lastLineStartInFile);
				file << indent << out << std::flush;
				//file << '\xd' << indent << out << std::flush;
			}
			if(isVerbose) std::cout << '\xd' << indent << out << std::flush;
    }

	template<typename T, typename U, typename V>
	void listOverFlush(T first, U middle, V last){
			if(isFile){
				file.seekp(lastLineStartInFile);
				file << indent << first << middle << last << std::flush;
				//file << '\xd' << indent << first << middle << last << std::flush;
			}
			if(isVerbose) std::cout << '\xd' << indent << first << middle << last << std::flush;
    }

	template<typename T, typename U, typename V, typename W, typename X>
	void listOverFlush(T first, U second, V third, W fourth, X fifth){
			if(isFile){
				file.seekp(lastLineStartInFile);
				file << indent << first << second << third << fourth << fifth << std::flush;
				//file << '\xd' << indent << first << second << third << fourth << fifth << std::flush;
			}
			if(isVerbose) std::cout << '\xd' << indent << first << second << third << fourth << fifth << std::flush;
    }

	//---------------------------------------------------------
	//fixed width
	//---------------------------------------------------------

	template<typename T>
	void flushFixedWidth(T out, int width){
		if(isFile) file << std::setw(width) << out << std::flush;
		if(isVerbose) std::cout << std::setw(width) << out << std::flush;
    }

	template<typename T>
	void flushNumberFixedWidth(T out, int precision, int width){
		if(isFile){
			std::ios::fmtflags old_settings = file.flags();
			file << std::setw(width) << std::fixed << std::setprecision(precision) << out << std::flush;
			file.flags(old_settings);
		}
		if(isVerbose){
			std::ios::fmtflags old_settings = std::cout.flags();
			std::cout << std::setw(width) << std::fixed << std::setprecision(precision) << out << std::flush;
			std::cout.flags(old_settings);
		}
    }

	template<typename T>
	void writeFixedWidth(T out, int width){
		if(isFile){
			file << std::setw(width) << out << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose) std::cout << std::setw(width) << out << std::endl;
    }

	template<typename T>
	void flushScientific(T out, int precision, int width){
		if(isFile){
			std::ios::fmtflags old_settings = file.flags();
			file << std::setw(width) << std::scientific << std::setprecision(precision)  << out << std::flush;
			file.flags(old_settings);
		}
		if(isVerbose){
			std::ios::fmtflags old_settings = std::cout.flags();
			std::cout << std::setw(width) << std::scientific << std::setprecision(precision)  << out << std::flush;
			std::cout.flags(old_settings);
		}
    }

    template<typename T>
	void writeScientific(T out, int precision){
		if(isFile){
			file << std::scientific << std::setprecision(precision) << out << std::endl;
			lastLineStartInFile=file.tellp();
		}
		if(isVerbose){
			std::cout.precision(precision);
			std::cout << std::scientific << std::setprecision(precision) << out << std::endl;
		}
    }

};





#endif /* TLOG_H_ */
