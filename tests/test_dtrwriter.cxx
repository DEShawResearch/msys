#include <stdio.h>
#include "molfile/dtrplugin.hxx"
#include "molfile/dtrframe.hxx"
#include "molfile/molfile_plugin.h"

#define FRAMES 10

int main(int argc, char *argv[]) {

    //
    // Write some frames to DTRs and ETRs
    //
    {
	//
	// This is a traditional DTR without a meta frame.
	//
	desres::molfile::DtrWriter dtr_no_meta(0, 100);
	dtr_no_meta.init("/tmp/test_no_meta.dtr", desres::molfile::DtrWriter::CLOBBER);
	
	//
	// This is a traditional DTR with a meta frame
	//
	desres::molfile::dtr::KeyMap meta_map;
	meta_map["PROVENANCE"] = desres::molfile::dtr::Key("test_etr", 9, desres::molfile::dtr::Key::TYPE_CHAR, false);
	desres::molfile::DtrWriter dtr_with_meta(0, 100);
	dtr_with_meta.init("/tmp/test_with_meta.dtr", desres::molfile::DtrWriter::CLOBBER, desres::molfile::DtrWriter::DTR, &meta_map);

	//
	// This is an ETR
	//
	desres::molfile::DtrWriter etr(0, 100);
	etr.init("/tmp/test_etr.etr", desres::molfile::DtrWriter::CLOBBER, desres::molfile::DtrWriter::ETR);
	
	//
	// Create some frames to add to each DTR/ETR
	//
	desres::molfile::dtr::KeyMap key_map;
	
	int single_int = 123;
	const char *string = "hello world";
	double chem_time;
	
	desres::molfile::dtr::Key single_int_key(&single_int, 1, desres::molfile::dtr::Key::TYPE_INT32, false);
	desres::molfile::dtr::Key string_key(string, 12, desres::molfile::dtr::Key::TYPE_CHAR, false);
	desres::molfile::dtr::Key chem_time_key(&chem_time, 1, desres::molfile::dtr::Key::TYPE_FLOAT64, false);
	
	for (double t = 0.0; t < (double) FRAMES; t += 1.0) {
	    chem_time = t;

	    key_map["SINGLE_INT"] = single_int_key;
	    key_map["STRING"] = string_key;
	    key_map["CHEM_TIME"] = chem_time_key;
	    
	    dtr_no_meta.append(t, key_map);
	    dtr_with_meta.append(t, key_map);
	    etr.append(t, key_map);
	}
    }

    {
	//
	// This is a traditional DTR without a meta frame.
	//
	desres::molfile::DtrWriter dtr_no_meta(0, 100);
	dtr_no_meta.init("/tmp/test_no_meta.dtr", desres::molfile::DtrWriter::APPEND);
	
	//
	// This is a traditional DTR with a meta frame
	//
	desres::molfile::dtr::KeyMap meta_map;
	meta_map["PROVENANCE"] = desres::molfile::dtr::Key("test_etr", 9, desres::molfile::dtr::Key::TYPE_CHAR, false);
	desres::molfile::DtrWriter dtr_with_meta(0, 100);
	dtr_with_meta.init("/tmp/test_with_meta.dtr", desres::molfile::DtrWriter::APPEND, desres::molfile::DtrWriter::DTR, &meta_map);

	//
	// This is an ETR
	//
	desres::molfile::DtrWriter etr(0, 100);
	etr.init("/tmp/test_etr.etr", desres::molfile::DtrWriter::APPEND, desres::molfile::DtrWriter::ETR);
	
	//
	// Create some frames to add to each DTR/ETR
	//
	desres::molfile::dtr::KeyMap key_map;
	
	int single_int = 123;
	const char *string = "hello world";
	double chem_time;
	
	desres::molfile::dtr::Key single_int_key(&single_int, 1, desres::molfile::dtr::Key::TYPE_INT32, false);
	desres::molfile::dtr::Key string_key(string, 12, desres::molfile::dtr::Key::TYPE_CHAR, false);
	desres::molfile::dtr::Key chem_time_key(&chem_time, 1, desres::molfile::dtr::Key::TYPE_FLOAT64, false);
	
	for (double t = FRAMES; t < (double) (FRAMES * 2.0); t += 1.0) {
	    chem_time = t;

	    key_map["SINGLE_INT"] = single_int_key;
	    key_map["STRING"] = string_key;
	    key_map["CHEM_TIME"] = chem_time_key;
	    
	    dtr_no_meta.append(t, key_map);
	    dtr_with_meta.append(t, key_map);
	    etr.append(t, key_map);
	}
    }

    //
    // Read back the frames
    //
    molfile_timestep_t ts;
    void *frame_data = NULL;
    desres::molfile::DtrReader etr("/tmp/test_etr.etr");
    etr.init();
    for (uint32_t frame_index = 0; frame_index < (2 * FRAMES); frame_index++) {
	double chem_time;
	desres::molfile::dtr::KeyMap keys = etr.frame(frame_index, &ts, &frame_data);
	keys["CHEM_TIME"].get(&chem_time);
	printf("CHEM_TIME from frame %u %f\n", frame_index, chem_time);
    }

    desres::molfile::DtrReader atr("/d/vault/envault-2/anton2/167/6925223/checkpoint.atr");
    atr.init();
    desres::molfile::dtr::KeyMap keys = atr.frame(0, NULL, &frame_data);
    char *cp = (char *) keys["FORMAT"].data;
    printf("atr FORMAT string %s\n", cp);
    
    return 0;
}
