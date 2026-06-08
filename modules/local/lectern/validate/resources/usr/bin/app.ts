#!/usr/bin/env node
//  Copyright (c) 2019, Ontario Institute for Cancer Research (OICR).

//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Affero General Public License as published
//  by the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Affero General Public License for more details.

//  You should have received a copy of the GNU Affero General Public License
//  along with this program. If not, see <https://www.gnu.org/licenses/>.

//  Author:Edmund Su <edmund.su@oicr.on.ca>
//     Linda Xiang <linda.xiang@oicr.on.ca>

import * as lectern from '@overture-stack/lectern-client';
import csvParser from 'csv-parser';
import fs from 'fs';
import path from 'path';
import minimist from 'minimist';


//node app.ts 
// --url https://dictionary-manager.submission.genomelibrary.ca 
// --directory . --study PCGLST0002


//Declare user defined variables
interface Args {
  path: string;
  url: string;
  dictionary: string;
  molecular : boolean;
}

const args = minimist<Args>(process.argv.slice(2));

const { url, dictionary, tsv , version, molecular} = args;
// If not part of molecular workflow, do not suppress participant foreign key check
const remove_participant = molecular ?? false;
// Read TSV files as an array
const tsvFiles: string[] = Array.isArray(tsv) ? tsv : tsv ? [tsv] : [];

//Error out if these keys are not supplied
if (!url || !dictionary || tsvFiles.length === 0) {
  console.error('Error: missing required arguments');
  console.error('Usage: ts app.ts --url https://dictionary-manager.submission.genomelibrary.ca --study \'PCTLST000000 custom schema\' --version 1.1 --tsv tsv1.tsv tsv2.tsv tsv3.tsv');
  process.exit(1);
}

//Read TSVs from filePath, removing null values, and save under entity named after file
const readTSV = (file: string): Promise<object[]> =>
  new Promise((resolve, reject) => {
    const records: object[] = [];
    fs.createReadStream(file)
      .pipe(csvParser({ separator: '\t' }))
      .on('data', (row) => {
        const filtered = Object.fromEntries(
          Object.entries(row)
            .filter(([_, v]) => v !== '' && v !== null && v !== undefined)
            .map(([k, v]) => {
              const num = Number(v);
              return [k, !isNaN(num) && v !== '' ? num : v];
            })
        );
        records.push(filtered);
      })
      .on('end', () => resolve(records))
      .on('error', reject);
  });

// Iterate through list and read TSVs saving to object name
const loadAllTSVs = async (files: string[]): Promise<Record<string, object[]>> => {
  const result: Record<string, object[]> = {};

  await Promise.all(
    files.map(async (file) => {
      const name = path.basename(file, '.tsv');
      result[name] = await readTSV(file);
    })
  );

  return result;
};

// Remove submitter_participant_id foreign key restriction
const removeParticipantRestrictions = (response: typeof apiResponse) => {
  return {
    ...response,
    data: {
      ...response.data,
      schemas: response.data.schemas.map((schema) => {
        const filteredForeignKeys = (schema.restrictions?.foreignKey ?? []).filter(
          (fk: any) =>
            fk.schema !== 'participant' &&
            !fk.mappings.some((m: any) => m.local === 'submitter_participant_id' || m.foreign === 'submitter_participant_id')
        );

        return {
          ...schema,
          restrictions: {
            ...schema.restrictions,
            foreignKey: filteredForeignKeys,
          },
        };
      }),
    },
  };
};

//Look for existence of dictionary
const filteredByNameDictionariesResult = await lectern.rest.listDictionaries(url, { name: dictionary });

//If dictionary does not exist, exit
if (filteredByNameDictionariesResult.data.length===0){
  console.error("Dictionary '"+dictionary+"' was not found");
  process.exit(1);
}

//Using latest version if version is not specified
let dictionaryVersion : string
if (!version){
  dictionaryVersion = filteredByNameDictionariesResult.data[filteredByNameDictionariesResult.data.length-1].version;
} else {
  dictionaryVersion = version
}


//Return dictionary schema
console.log("Using "+dictionary+" "+ dictionaryVersion)
const pcglDictionary = await lectern.rest.getDictionary(url, {
    name: dictionary,
    version: dictionaryVersion,
});

//Load TSV data
const data = await loadAllTSVs(tsvFiles);

// If molecular suppress submitter participant ID restriction
let updated_dictionary: unknown;

if (remove_participant) {
  updated_dictionary = removeParticipantRestrictions(pcglDictionary);
} else {
  updated_dictionary = pcglDictionary;
}

//Validate data through schema
const processingResult = lectern.validate.validateDictionary(data, updated_dictionary.data);

//Return errors if detected
if (processingResult.valid){
    process.exit(0)
} else {
    for (const detail in processingResult.details){
        for (const record in processingResult.details[detail].invalidRecords){
            for (const errorRecord in processingResult.details[detail].invalidRecords[record].recordErrors){
              console.error(
                  JSON.stringify(
                      {
                          index: processingResult.details[detail].invalidRecords[record].recordIndex,
                          errors: processingResult.details[detail].invalidRecords[record].recordErrors[errorRecord].errors,
                          reason: processingResult.details[detail].invalidRecords[record].recordErrors[errorRecord].reason,
                          fieldName: processingResult.details[detail].invalidRecords[record].recordErrors[errorRecord].fieldName,
                          fieldValue: processingResult.details[detail].invalidRecords[record].recordErrors[errorRecord].fieldValue,
                      },
                      null,
                      2
                  )
              );
            }
        }
    }
    console.log("Validation Failed")
    process.exit(1);
}
