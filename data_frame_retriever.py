from collections import defaultdict
from Core.store.postgres_query_manager import PostgresQueryManager
from Insights.util.get_all_insights_info import get_all_insights_info
import pandas as pd

TAGS_QUERY = """
    with recent_articles as (
        SELECT
            id,
            process_date
        FROM
            tbl_article_metadata
        WHERE
            process_date >= now() - %(time_span)s::interval
    ),
    recent_tags as (
        SELECT
            company_id as resource_id,
            recent_articles.id as article_id,
            CAST(recent_articles.process_date as date) as article_process_day,
            CAST(recent_articles.process_date - '1 week'::interval as date) as week_ago
        FROM
            tbl_article_company_tags
        JOIN
            tbl_article_company_tag_mentions
        ON tbl_article_company_tag_mentions.company_tag_id = tbl_article_company_tags.id
        JOIN 
            recent_articles
        ON recent_articles.id = tbl_article_company_tags.article_id
        UNION ALL
        SELECT
            person_id as resource_id,
            recent_articles.id as article_id,
            CAST(recent_articles.process_date as date) as article_process_day,
            CAST(recent_articles.process_date - '1 week'::interval as date) as week_ago
        FROM
            tbl_article_person_tags
        JOIN
            tbl_article_person_tag_mentions
        ON tbl_article_person_tag_mentions.person_tag_id = tbl_article_person_tags.id
        JOIN 
            recent_articles
        ON recent_articles.id = tbl_article_person_tags.article_id
    )
    SELECT
        tbl_primitive_label_resources.label_id,
        recent_tags.article_id,
        recent_tags.article_process_day as date,
        CAST(recent_tags.week_ago AS TEXT)
    FROM
        recent_tags
    JOIN
        tbl_primitive_label_resources
    ON tbl_primitive_label_resources.resource_id = recent_tags.resource_id
"""

UNMATCHED_MENTIONS_QUERY = """
    with recent_articles as (
        SELECT
            id,
            process_date
        FROM
            tbl_article_metadata
        WHERE
            process_date >= now() - %(time_span)s::interval
    )
    SELECT
        -1 as label_id,
        text as mention_text,
        article_id,
        CAST(recent_articles.process_date as date) as date,
        CAST(recent_articles.process_date - '1 week'::interval as date) as week_ago,
        mention_type as ner_type
    FROM
        tbl_article_unmatched_mention_tags
    JOIN
        recent_articles
    ON recent_articles.id = tbl_article_unmatched_mention_tags.article_id
"""

VW_LABELS_QUERY = """
    SELECT 
        COALESCE(vw_labels.arkera_name,vw_labels.wikidata_name,vw_labels.google_name) as label_name, 
        id as label_id,
        resource_id,
        resource_type
    FROM vw_labels 
    WHERE resource_type in ('companies','persons','organizations')
"""


class InsightsMonitoringDataFrameRetriever:

    def __init__(self, queryManager):
        self.controlSheetsDF = get_all_insights_info().get_entities_of_interest()[['label_id', 'insightsTopic', 'insightsCountry']]
        self.queryManager = queryManager

    def get_data_frame(self, days=9):
        insightsTagsDF = pd.merge(self.controlSheetsDF, self.queryManager.run_query(query=TAGS_QUERY, params={'time_span': f'{days} days'}), on='label_id', how='left').replace({pd.np.nan: 'N/A'})
        unlinkedMentionsDF = self.queryManager.run_query(query=UNMATCHED_MENTIONS_QUERY, params={'time_span': f'{days} days'})
        insightsTagsDF = insightsTagsDF.append(unlinkedMentionsDF).replace({pd.np.nan: 'N/A'})
        groupByColumns = [x for x in insightsTagsDF.columns.tolist() if x not in {'article_id'}]
        finalDF = insightsTagsDF.groupby(groupByColumns).apply(lambda df: pd.Series({
            'count': len(df.query('article_id != "N/A"')),
            'article_ids': [int(x) for x in df.query('article_id != "N/A"')['article_id']]
        })).reset_index()
        mergedControlSheetsDF = pd.merge(finalDF, self.queryManager.run_query(query=VW_LABELS_QUERY), how='left', on='label_id')
        return mergedControlSheetsDF.sort_values(['count'], ascending=False)
