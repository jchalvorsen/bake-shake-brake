#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include "boost/iostreams/filter/zlib.hpp"


#define READ_BYTES3(var, buf) *((char*)&var+2) = (buf)[0]; *((char*)&var+1) = (buf)[1]; *((char*)&var+0) = (buf)[2];
#define READ_BYTES4(var, buf) *((char*)&var+3) = (buf)[0]; *((char*)&var+2) = (buf)[1]; *((char*)&var+1) = (buf)[2]; *((char*)&var+0) = (buf)[3];
#define IS_STAIRS(block) ((block)==53 || (block)==67 || (block)==108 || (block)==109 || (block)==114)
#define IS_SLAB(block) ((block)==44)
#define IS_DECORATION(block) ((block)==50 || (block)==31 || (block)==72 || ((block)>36 && (block)<41) || (block)==55 || (block)==63)
#define IS_DOOR(block) ((block)==64 || (block)==71)

using namespace std;
using namespace boost::iostreams;

int main(int argc, char **argv) {

	// request bounding box
	int xmin, xmax;
	int ymin, ymax;
	int zmin, zmax;
	string regionFolder;
	if(argc>1) {
		xmin = atoi(argv[1]);
		xmax = atoi(argv[2]);
		ymin = atoi(argv[3]);
		ymax = atoi(argv[4]);
		zmin = atoi(argv[5]);
		zmax = atoi(argv[6]);
		if(argc > 7)
			regionFolder = argv[7];
		else
			regionFolder = "region";
	} else {
		cout << "xmin: "; cin >> xmin;
		cout << "xmax: "; cin >> xmax;
		cout << "ymin: "; cin >> ymin;
		cout << "ymax: "; cin >> ymax;
		cout << "zmin: "; cin >> zmin;
		cout << "zmax: "; cin >> zmax;
		regionFolder = "region";
	}

	// allocate result data arrays
	int dx = xmax-xmin;
	int dy = ymax-ymin;
	int dz = zmax-zmin;
	unsigned char ***minecraft_block = new unsigned char**[dx];
	unsigned char ***minecraft_data  = new unsigned char**[dx]; // meta-information like stair orientation
	for(int i=0; i<dx; i++) {
		minecraft_block[i] = new unsigned char*[dy];
		minecraft_data[i]  = new unsigned char*[dy];
		for(int j=0; j<dy; j++) {
			minecraft_block[i][j] = new unsigned char[dz];
			minecraft_data[i][j]  = new unsigned char[dz];
		}
	}

	/***************************************************************************************************************
	
	*******************         Read the raw data from region .mca-files                        ********************

	***************************************************************************************************************/

	// find chunk span and region file(s)
	int chunkX_min  = (int) floor((double)  xmin   / 16);
	int chunkX_max  = (int) floor((double)(xmax-1) / 16);
	int chunkZ_min  = (int) floor((double)  zmin   / 16);
	int chunkZ_max  = (int) floor((double)(zmax-1) / 16);
	int regionX_min = (int) floor((double)chunkX_min / 32);
	int regionZ_min = (int) floor((double)chunkZ_min / 32);
	int regionX_max = (int) floor((double)chunkX_max / 32);
	int regionZ_max = (int) floor((double)chunkZ_max / 32);

	char dataBuffer[0x800];
	char blockBuffer[0x1000];

	// loop over all region files
	for(int regionX=regionX_min; regionX<=regionX_max; regionX++) {
		for(int regionZ=regionZ_min; regionZ<=regionZ_max; regionZ++) {

			int loc_chunkX_min = (chunkX_min <    regionX  *32) ?    0  : chunkX_min % 32;
			int loc_chunkX_max = (chunkX_max >= (regionX+1)*32) ?   31  : chunkX_max % 32;
			int loc_chunkZ_min = (chunkZ_min <    regionZ  *32) ?    0  : chunkZ_min % 32;
			int loc_chunkZ_max = (chunkZ_max >= (regionZ+1)*32) ?   31  : chunkZ_max % 32;
			// some architechtures give negative modulus results
			if(loc_chunkX_min<0) loc_chunkX_min+=32;
			if(loc_chunkX_max<0) loc_chunkX_max+=32;
			if(loc_chunkZ_min<0) loc_chunkZ_min+=32;
			if(loc_chunkZ_max<0) loc_chunkZ_max+=32;
			
			stringstream buf;
			buf << regionFolder;
			buf << "/r." << regionX << "." << regionZ << ".mca";
			string filename = buf.str();
			cout << "Parsing region file " << filename << endl;

			// loop over all chunks in this region file
			for(int chunkX=loc_chunkX_min; chunkX<=loc_chunkX_max; chunkX++) {
				for(int chunkZ=loc_chunkZ_min; chunkZ<=loc_chunkZ_max; chunkZ++) {

					// open the region file
					ifstream file (filename.c_str(), ios::in | ios::binary);
					if( !file.is_open() ) {
						cerr << "Error opening file " << filename << endl;
						exit(1);
					}
					cout << "  Reading file " << endl;

					cout << "  Reading chunk (" << chunkX << ", " << chunkZ << ")\n";

					// read the chunk header
					int offsetX = (chunkX<0) ? chunkX%32+32 : chunkX%32; // some silly architecthures give negative modulus results
					int offsetZ = (chunkZ<0) ? chunkZ%32+32 : chunkZ%32;
					int chunkheader_offset = 4*(offsetX + offsetZ * 32);
					file.seekg(chunkheader_offset);
					char buffer[256];
					unsigned int chunk_position = 0;
					file.read(buffer, 4);
					READ_BYTES3(chunk_position, buffer);
					if(chunk_position < 2) {
						cerr << "Error: Data not generated for current chunk (check proper bounding box)\n";
						file.close();
						exit(2);
					}

					// read the chunk encrypted block
					file.seekg(chunk_position*0x1000);
					unsigned int binary_size;
					unsigned int compression_type;
					file.read(buffer, 5);
					READ_BYTES4(binary_size, buffer);
					compression_type = buffer[4];
					if(compression_type != 2) {
						cerr << "Error reading map file. Unsupported encryption type: " << compression_type << "\n";
						file.close();
						exit(3);
					}

					// decrypt the encrypted block 
					file.seekg(chunk_position*0x1000+5);
					filtering_streambuf<input> in;
					in.push(zlib_decompressor()); 
					in.push(file);

					// debug code, dumping data to files
					// cout << "Writing the decompressed chunk (" << chunkX << ", " << chunkZ << ") to oneChunk.out" << endl;
					// char outputFileName[256];
					// sprintf(outputFileName, "oneChunk-%d-%d.out", chunkX, chunkZ);
					const char *outputFileName = "oneChunk.out";
					ofstream outfile;
					outfile.open(outputFileName, ios::binary);

					// write decrypted stream to new file
					boost::iostreams::copy(in, outfile);
					outfile.close();

					// dump compressed data to file
					// cout << "Writing the same compressed chunk in compressed.gz of " << endl;
					// ofstream compressedFile;
					// compressedFile.open("compressed.gz", ios::binary);
					// file.seekg(chunk_position*0x1000+5);
					// boost::iostreams::copy(file, compressedFile);
					// compressedFile.close();

					// read the chunk from file
					ifstream is("oneChunk.out", ifstream::binary);

					bool found = false;
					char oneChar = get(is);
					while(!found && !is.eof()) {
						// scan chunk for 'Sections' keyword
						while(oneChar != 'S' && !is.eof())
							oneChar = get(is);
						read(is, buffer, 7);
						buffer[7] = 0;
						if(strncmp(buffer, "ections", 7) == 0) {
							found = true;
							// cout << "Sections found!!!\n";
						} else {
							putback(is, 7);
							oneChar = get(is);
						}
					}
					
					// loop over all Y-layers (max 16)
					int chunkY = -1;
					while(!is.eof()) {
						
						fill(dataBuffer,  dataBuffer+0x800,   0);
						fill(blockBuffer, blockBuffer+0x1000, 0);
					
						// scan chunk for 'Data' keyword
						found = false;
						while(!found && !is.eof()) {
							while(oneChar != 'D' && !is.eof())
								oneChar = get(is);
							read(is, buffer, 3);
							buffer[3] = 0;
							if(strncmp(buffer, "ata", 3) == 0) {
								found = true;
								// cout << "Chunk meta data sucessfully detected\n";

								// skip some header information 
								read(is, buffer, 4);
								// read main block of data
								read(is, dataBuffer, 0x800);
							} else {
								putback(is, 3);
								oneChar = get(is);
							}
						}

						found = false;
						while(!found && !is.eof()) {
							// scan chunk for 'Y....Blocks' keyword
							while(oneChar != 'Y' && !is.eof())
								oneChar = get(is);
							read(is, buffer, 10);
							buffer[11] = 0;
							if(strncmp(buffer+4, "Blocks", 6) == 0) {
								chunkY = buffer[0];
								found = true;

								// cout << "Chunk block data sucessfully detected\n";

								// skip some header information 
								read(is, buffer, 4);
								// read main block of data
								read(is, blockBuffer, 0x1000);

							} else {
								putback(is, 10);
								oneChar = get(is);
							}
						}

						if(is.eof())
							continue;

						cout << "    Y layer : " << chunkY << endl;
						if( ymin >= (chunkY+1)*16 || ymax < chunkY*16 )
							continue;

						int xSize = 16;
						int ySize = 16;
						int zSize = 16;
						int globXmin = regionX*32*16 + chunkX*16;
						int globYmin =                 chunkY*16;
						int globZmin = regionZ*32*16 + chunkZ*16;
						int locXmin = (xmin <   globXmin   )  ?    0  : xmin % 16;
						int locXmax = (xmax >=  globXmin+16 ) ? xSize : xmax % 16;
						int locYmin = (ymin <   globYmin   )  ?    0  : ymin % 16;
						int locYmax = (ymax >=  globYmin+16 ) ? ySize : ymax % 16;
						int locZmin = (zmin <   globZmin   )  ?    0  : zmin % 16;
						int locZmax = (zmax >=  globZmin+16 ) ? zSize : zmax % 16;

						// some architechtures give negative modulus results
						if(locXmin<0) locXmin+=16;
						if(locXmax<0) locXmax+=16;
						if(locYmin<0) locYmin+=16;
						if(locYmax<0) locYmax+=16;
						if(locZmin<0) locZmin+=16;
						if(locZmax<0)locZmax+=16;
						unsigned char oneBlock;
						int k=0;
						for(int y=0; y<ySize; y++) {
							for(int z=0; z<zSize; z++) {
								for(int x=0; x<xSize; x++) {
									oneBlock = blockBuffer[k++];
									if( locXmin <= x && x < locXmax &&
										locYmin <= y && y < locYmax &&
										locZmin <= z && z < locZmax )
										minecraft_block[globXmin+x-xmin][globYmin+y-ymin][globZmin+z-zmin] = oneBlock;
								}
							}
						}
						// cout << "Chunk block ID sucessfully read\n";

						unsigned char oneByte;
						k=0;
						for(int y=0; y<ySize; y++) {
							for(int z=0; z<zSize; z++) {
								for(int x=0; x<xSize; x+=2) {
									oneByte = dataBuffer[k++];
									if( locXmin <= x && x < locXmax &&
										locYmin <= y && y < locYmax && 
										locZmin <= z && z < locZmax) {
										if(xmin <= globXmin+x && globXmin+x < xmax ) 
											minecraft_data[globXmin+x-xmin  ][globYmin+y-ymin][globZmin+z-zmin] = (oneByte & 0xF); // least significant 4 bits first
										if(xmin <= globXmin+x+1 && globXmin+x+1 < xmax ) 
											minecraft_data[globXmin+x-xmin+1][globYmin+y-ymin][globZmin+z-zmin] = (oneByte >> 4);
										/*
										if(oneByte != 0)  {
												printf("nonzero data byte read: %3u\n", oneByte);
												printf("  x      = %d\n  y      = %d\n  z      = %d\n",
												globXmin+x-xmin,
												globYmin+y-ymin,
												globZmin+z-zmin);
											if(xmin <= globXmin+x && globXmin+x < xmax ) 
												printf("  value1 = %u (block ID: %3u)\n", minecraft_data[ globXmin+x-xmin  ][globYmin+y-ymin][globZmin+z-zmin],
												                                          minecraft_block[globXmin+x-xmin  ][globYmin+y-ymin][globZmin+z-zmin]);
											if(xmin <= globXmin+x+1 && globXmin+x+1 < xmax ) 
												printf("  value2 = %u (block ID: %3u)\n", minecraft_data[ globXmin+x-xmin+1][globYmin+y-ymin][globZmin+z-zmin],
												                                          minecraft_block[globXmin+x-xmin+1][globYmin+y-ymin][globZmin+z-zmin]);
										}
										*/
									}
								}
							}
						}
						// cout << "Chunk meta data sucessfully read\n";
						
					}
					cout<< "  Done with chunk (" << chunkX << ", " << chunkZ << ")\n";

					// close file and make ready for a new one
					file.close();
				}
			}

		}
	}

	// verify by dumping data to screen
	/*
	int y ;
	cout << "y-plane: ";
	cin >> y;
	bool meta = y<0;
	if(y<0)
		y = -y;
	while(ymin <= y && y < ymax) {
		for(int z=dz; z-->0; ) {
			for(int x=dx; x-->0; ) {
				if(meta)
					printf("%3u ", minecraft_data[x][y-ymin][z]);
				else
					printf("%3u ", minecraft_block[x][y-ymin][z]);
			}
			cout << endl;
		}
		cout << "y-plane: ";
		cin >> y;
		meta=y<0;
		if(y<0)
			y = -y;
	}
	*/


	/***************************************************************************************************************
	
	*******************         Build the finite element data structures                        ********************

	***************************************************************************************************************/

	// allocate memory blocks
	int ****finiteElement = new int***[2*dx];
	for(int i=0; i<2*dx; i++) {
		finiteElement[i] = new int**[2*dy];
		for(int j=0; j<2*dy; j++) {
			finiteElement[i][j] = new int*[2*dz];
			for(int k=0; k<2*dz; k++) {
				finiteElement[i][j][k] = new int[11];
				for(int l=0; l<11; l++) // 11 meta-values. The global index, 8 corner nodes, one for the material and one for meta-data
					finiteElement[i][j][k][l] = -1;
			}
		}
	}
	double ****node = new double***[2*dx+1];
	for(int i=0; i<2*dx+1; i++) {
		node[i] = new double**[2*dy+1];
		for(int j=0; j<2*dy+1; j++) {
			node[i][j] = new double*[2*dz+1];
			for(int k=0; k<2*dz+1; k++) {
				node[i][j][k] = new double[4]; // 4 meta-values, the global index and the (x,y,z)-coordinates 
				for(int l=0; l<4; l++)
					node[i][j][k][l] = -1;
			}
		}
	}
	int nodeCount    = 0;
	int elementCount = 0;

	for(int i=0; i<2*dx; i++) {
		for(int j=0; j<2*dy; j++) {
			for(int k=0; k<2*dz; k++) {
				if(minecraft_block[i/2][j/2][k/2] != 0 &&
				   minecraft_block[i/2][j/2][k/2] != 9) {
					if(IS_DECORATION(minecraft_block[i/2][j/2][k/2]))         // skip decorations altogether
						continue;
					if(IS_SLAB(minecraft_block[i/2][j/2][k/2]) && j%2==1)     // skip upper half of slabs
						continue;
					if(IS_STAIRS(minecraft_block[i/2][j/2][k/2]) ) { // skip the proper part of stair blocks
						unsigned char direction   = minecraft_data[i/2][j/2][k/2] & 0x3;
						unsigned char orientation = (minecraft_data[i/2][j/2][k/2] >> 2) & 0x1; // 4th bit 1 = upside-down stairs
						if(direction==0 && i%2==0 && j%2!=orientation) // positive x
							continue;
						if(direction==1 && i%2==1 && j%2!=orientation) // negative x
							continue;
						if(direction==2 && k%2==0 && j%2!=orientation) // positive z
							continue;
						if(direction==3 && k%2==1 && j%2!=orientation) // negative z
							continue;
					}
					if(IS_DOOR(minecraft_block[i/2][j/2][k/2])) { // skip vertical half of doors
						unsigned char info = minecraft_data[i/2][j/2][k/2];
						bool swung = info & 0x4;
						info &= 0x3; // clear swung & top information
						if((info == 0 && !swung) || (info == 3 && swung)) {// west position
							if(i%2==1) continue;
						} else if((info == 1 && !swung) || (info == 0 && swung)) { // north position
							if(k%2==1) continue;
						} else if((info == 2 && !swung) || (info == 1 && swung)) { // east position
							if(i%2==0) continue;
						} else if((info == 3 && !swung) || (info == 2 && swung)) { // south position
							if(k%2==0) continue;
						}
					}
					finiteElement[i][j][k][0] = elementCount++;

					// loop over 8 corner points of this element
					int cornerCount = 1;
					for(int w=0; w<2; w++) {
						for(int v=0; v<2; v++) {
							for(int u=v; 0<=u && u<2; (v)?u--:u++) { // loops (u,v,w) indices by (0,0,0), (0,1,0), (1,1,0), (1,0,0) ...
								if(node[i+u][j+v][k+w][0] == -1) {
									node[i+u][j+v][k+w][0] = nodeCount++; 
									node[i+u][j+v][k+w][1] = ((double) i+u )/2.0 + xmin;
									node[i+u][j+v][k+w][2] = ((double) j+v )/2.0 + ymin;
									node[i+u][j+v][k+w][3] = ((double) k+w )/2.0 + zmin;
								}
								finiteElement[i][j][k][cornerCount++] = node[i+u][j+v][k+w][0];
							}
						}
					}
					// set the material- and meta-properties of this element
					finiteElement[i][j][k][9]  = minecraft_block[i/2][j/2][k/2];
					finiteElement[i][j][k][10] = (int) minecraft_data[i/2][j/2][k/2];
				}

/*
				if(minecraft_data[i/2][j/2][k/2] != 0) {
					cout << "Meta data != 0 for (" << i << ", " << j << ", " << k << ")\n";
					cout << "  Block: " << minecraft_block[i/2][j/2][k/2] << endl;
					cout << "  Meta:  " << minecraft_data[i/2][j/2][k/2] << endl;
				}
*/

			}
		}
	}

	// build an ordered set of the nodes (they're scattered in the node-array wrt global indexing)
	int **nodeIndex = new int*[nodeCount];
	for(int i=0; i<2*dx+1; i++) {
		for(int j=0; j<2*dy+1; j++) {
			for(int k=0; k<2*dz+1; k++) {
				int globi = node[i][j][k][0];
				if(globi != -1) {
					nodeIndex[globi] = new int[3];
					nodeIndex[globi][0] = i;
					nodeIndex[globi][1] = j;
					nodeIndex[globi][2] = k;
				}
			}
		}
	}
	
	/***************************************************************************************************************
	
	*******************         Dump results to files                                           ********************

	***************************************************************************************************************/

	ofstream vtf_file;
	vtf_file.open("minecraft_blocks.vtf");
	if(!vtf_file.good()) {
		cerr << "Error: opening \"minecraft_blocks.vtf\"";
		exit(4);
	}
	cout << "writing file \"minecraft_blocks.vtf\"\n";

	vtf_file << "*VTF-1.00\n";
	vtf_file << "\n";
	vtf_file << "*INTERNALSTRING 40001\n";
	vtf_file << "VTF Writer Version info:\n";
	vtf_file << "APP_INFO: GLview Express Writer: 1.1-12\n";
	vtf_file << "GLVIEW_API_VER: 2.1-22\n";
	vtf_file << "EXPORT_DATE: 2011-10-19 12:48:00\n";
	vtf_file << "\n";
	vtf_file << "*NODES 1 \n";
	for(int n=0; n<nodeCount; n++) {
		int i = nodeIndex[n][0];
		int j = nodeIndex[n][1];
		int k = nodeIndex[n][2];
		vtf_file << "  " <<  node[i][j][k][1] << " " << node[i][j][k][2] << " " << node[i][j][k][3]  << endl;
	}
	vtf_file << "\n";
	vtf_file << "*ELEMENTS 1\n";
	vtf_file << "%NODES #1\n";
	vtf_file << "%NAME \"Patch 1\"\n";
	vtf_file << "%NO_ID\n";
	vtf_file << "%MAP_NODE_INDICES\n";
	vtf_file << "%PART_ID 1\n";
	vtf_file << "%HEXAHEDRONS\n";

	for(int i=0; i<2*dx; i++) {
		for(int j=0; j<2*dy; j++) {
			for(int k=0; k<2*dz; k++) {
				if(finiteElement[i][j][k][0] != -1) {
					for(int l=1; l<9; l++) 
						vtf_file << (finiteElement[i][j][k][l]+1) << " "; // vtf-file is 1-indexing things
					vtf_file << endl;
				}
			}
		}
	}

	vtf_file << "\n";
	vtf_file << "*RESULTS 2\n";
	vtf_file << "%NO_ID\n";
	vtf_file << "%DIMENSION 1\n";
	vtf_file << "%PER_ELEMENT #1\n";

	for(int i=0; i<2*dx; i++)
		for(int j=0; j<2*dy; j++)
			for(int k=0; k<2*dz; k++)
				if(finiteElement[i][j][k][0] != -1)
					vtf_file << (finiteElement[i][j][k][9]) << endl; // dump block code

	vtf_file << "\n";
	vtf_file << "*RESULTS 3\n";
	vtf_file << "%NO_ID\n";
	vtf_file << "%DIMENSION 1\n";
	vtf_file << "%PER_ELEMENT #1\n";

	for(int i=0; i<2*dx; i++)
		for(int j=0; j<2*dy; j++)
			for(int k=0; k<2*dz; k++)
				if(finiteElement[i][j][k][0] != -1)
					vtf_file << (finiteElement[i][j][k][10]) << endl; // dump block code

	vtf_file << "\n";
	vtf_file << "*GLVIEWGEOMETRY 1\n";
	vtf_file << "%STEP 1\n";
	vtf_file << "%ELEMENTS\n";
	vtf_file << "1\n";
	vtf_file << "\n";
	vtf_file << "*GLVIEWSCALAR 1\n";
	vtf_file << "%NAME \"Block code\"\n";
	vtf_file << "%STEP 1\n";
	vtf_file << "2\n";
	vtf_file << "\n";
	vtf_file << "*GLVIEWSCALAR 2\n";
	vtf_file << "%NAME \"Block meta\"\n";
	vtf_file << "%STEP 1\n";
	vtf_file << "3\n";

	vtf_file.close();

	ofstream nodefile;
	nodefile.open("nodes.m");
	if(!nodefile.good()) {
		cerr << "Error: opening \"nodes.m\"";
		exit(4);
	}
	cout << "writing file \"nodes.m\"\n";
	for(int n=0; n<nodeCount; n++) {
		int i = nodeIndex[n][0];
		int j = nodeIndex[n][1];
		int k = nodeIndex[n][2];
		nodefile << node[i][j][k][1] << " " << node[i][j][k][2] << " " << node[i][j][k][3]  << endl;
	}
	nodefile.close();


	ofstream elementfile;
	elementfile.open("elements.m");
	if(!elementfile.good()) {
		cerr << "Error: opening \"elements.m\"";
		exit(4);
	}
	cout << "writing file \"elements.m\"\n";
	for(int i=0; i<2*dx; i++) {
		for(int j=0; j<2*dy; j++) {
			for(int k=0; k<2*dz; k++) {
				if(finiteElement[i][j][k][0] != -1) {
					for(int l=1; l<11; l++) 
						elementfile << (finiteElement[i][j][k][l]) << " "; 
					elementfile << endl;
				}
			}
		}
	}
	elementfile.close();
	
}
